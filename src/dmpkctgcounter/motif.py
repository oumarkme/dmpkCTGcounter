"""Repeat motif analysis: purity, interruption blocks, and periodic runs.

The motif pass never changes a repeat count; it only qualifies reads (purity
filter), reports candidate interruptions, and measures how far a periodic
run extends into a read that has a single anchor.

Orientation: tract segments passed to ``analyze_tract`` are given on the
reference plus strand ((CAG)n) and are reverse-complemented here so that the
unit is CTG, the gene strand.
"""

from __future__ import annotations

import dataclasses
import re

from dmpkctgcounter import locus

UNIT = locus.GENE_MOTIF  # "CTG"
# Single-substitution variants of CTG and CAG (gene strand). Tokens that
# are not CTG are interruption candidates.
INTERRUPTION_TOKENS = frozenset({"CCG", "CGG", "CTC", "CAG", "TTG"})
_ALL_TOKENS = INTERRUPTION_TOKENS | {UNIT}
# Two same-type interruption tokens closer than this many unit positions
# form one block.
_BLOCK_MAX_SPACING = 4
RUN_WINDOW = 30
RUN_MIN_PURITY = 0.8
# Plus-strand forms of CTG and its one-base variants (reverse complements of
# the gene-strand tokens above): CAG, CGG, CCG, GAG, CTG, CAA.
_PLUS_TOKENS = frozenset(
    locus.reverse_complement(t) for t in INTERRUPTION_TOKENS | {UNIT}
)
_PLUS_RE = re.compile("|".join(sorted(_PLUS_TOKENS)))


@dataclasses.dataclass(frozen=True)
class Block:
    """A block of >= 2 same-type interruption tokens in one read.

    Attributes:
      token: Interruption token (gene strand).
      offset: Distance in unit positions from the 3' end of the tract (the
        gene-strand 3' anchor, the plus-strand left anchor) to the block's
        3'-most token. An absolute offset, so it does not depend on the
        length of the allele.
      count: Number of tokens in the block.
    """

    token: str
    offset: int
    count: int


@dataclasses.dataclass(frozen=True)
class TractAnalysis:
    """Result of tokenizing one inter-anchor segment.

    Attributes:
      purity: ``3 * n_CTG / length``, the strict CTG fraction.
      periodic: ``3 * (n_CTG + n_interruption tokens) / length``. Tokens
        differ from CTG by one base; this fraction measures how much of the
        segment is unit-like and is the one used to reject non-repeat reads.
      n_units: Number of CTG tokens.
      ttg_fraction: TTG tokens as a fraction of all unit-like tokens.
      blocks: Interruption blocks.
    """

    purity: float
    periodic: float
    n_units: int
    ttg_fraction: float
    blocks: tuple[Block, ...]


def analyze_tract(plus_segment: str) -> TractAnalysis:
    """Tokenizes a tract segment into CTG and interruption tokens.

    Scans left to right in steps of one triplet. A triplet that is neither
    CTG nor a known interruption is resynchronized by looking for a token
    at offset +1 or +2 (single insertions or deletions); otherwise it is
    skipped. Purity is ``3 * n_CTG / length``.

    Args:
      plus_segment: Bases between the two margins, plus strand (length is
        ``gap - 7`` in the read).

    Returns:
      Purity, number of CTG tokens, TTG fraction, and interruption blocks.
    """
    seg = locus.reverse_complement(plus_segment)
    n_total = len(seg)
    i = 0
    n_units = 0
    index = 0
    found: dict[str, list[int]] = {}
    while i + 3 <= n_total:
        tri = seg[i : i + 3]
        if tri == UNIT:
            n_units += 1
            index += 1
            i += 3
        elif tri in INTERRUPTION_TOKENS:
            found.setdefault(tri, []).append(index)
            index += 1
            i += 3
        else:
            for offset in (1, 2):
                if seg[i + offset : i + offset + 3] in _ALL_TOKENS:
                    i += offset
                    break
            else:
                index += 1
                i += 3
    purity = 3 * n_units / n_total if n_total else 0.0
    n_variant = sum(len(v) for v in found.values())
    periodic = 3 * (n_units + n_variant) / n_total if n_total else 0.0
    n_tokens = n_units + n_variant
    ttg_fraction = len(found.get("TTG", [])) / n_tokens if n_tokens else 0.0
    blocks = []
    for token, positions in found.items():
        chain = [positions[0]]
        for pos in [*positions[1:], None]:
            if pos is not None and pos - chain[-1] <= _BLOCK_MAX_SPACING:
                chain.append(pos)
                continue
            if len(chain) >= 2:
                offset3 = max(index - 1 - chain[-1], 0)
                blocks.append(Block(token, offset3, len(chain)))
            chain = [pos] if pos is not None else []
    return TractAnalysis(purity, periodic, n_units, ttg_fraction, tuple(blocks))


def periodic_run(seq: str, origin: int, direction: int) -> tuple[int, bool]:
    """Measures the unit-periodic run from ``origin`` into a read.

    Windows of 30 bp pass while at least 80% of their length is covered by
    non-overlapping unit-like tokens (CAG or a one-base variant on the plus
    strand; this tolerates single indels and the CTG to TTG changes seen in
    expanded repeats). After the last passing window the run is extended
    triplet by triplet.

    Args:
      seq: Read, plus strand (reference-forward).
      origin: 0-based offset where the repeat starts (after the margin).
      direction: +1 downstream (left anchor found), -1 upstream.

    Returns:
      Bases covered by the run, and whether the run reaches within one
      window of the read end (the read stops inside the repeat).
    """
    min_tokens = RUN_MIN_PURITY * RUN_WINDOW / 3
    pos = origin
    while True:
        if direction > 0:
            window = seq[pos : pos + RUN_WINDOW]
        else:
            window = seq[max(0, pos - RUN_WINDOW) : pos]
        if (
            len(window) < RUN_WINDOW
            or len(_PLUS_RE.findall(window)) < min_tokens
        ):
            break
        pos += direction * RUN_WINDOW
    while True:
        triplet = (
            seq[pos : pos + 3] if direction > 0 else seq[max(0, pos - 3) : pos]
        )
        if triplet not in _PLUS_TOKENS:
            break
        pos += direction * 3
    covered = abs(pos - origin)
    remaining = len(seq) - pos if direction > 0 else pos
    return covered, remaining < RUN_WINDOW
