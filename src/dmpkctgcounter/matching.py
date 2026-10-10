"""Anchor search and per-read repeat measurement.

Offsets are 0-based half-open positions in the read string. The CIGAR of a
record is never used for the repeat length: both anchors are located in the
read's own sequence (including soft-clipped bases) with edlib HW alignment,
and the repeat count is ``(gap - 7) / 3`` where ``gap`` is the distance
between the inner edges of the anchors.
"""

from __future__ import annotations

import dataclasses

import edlib

from dmpkctgcounter import locus, motif


@dataclasses.dataclass(frozen=True)
class Hit:
    """One anchor occurrence: ``read[start:end]`` with ``edits`` edits."""

    start: int
    end: int
    edits: int


@dataclasses.dataclass(frozen=True)
class Measure:
    """Result of measuring one read.

    Attributes:
      status: "spanning", "one_anchor" (periodic run reaches the read end),
        "partial" (one anchor, run does not reach the end), "rejected", or
        "no_anchor".
      reason: Reason code for rejected, partial, and no_anchor reads.
      gap: Bases between the inner anchor edges (spanning and rejected).
      units: Repeat units ``(gap - 7) / 3`` as a real number.
      remainder: ``(gap - 7) mod 3``, recorded and never used as a filter.
      purity: Strict CTG token purity of the tract (spanning reads).
      periodic: Fraction of the tract in unit-like tokens (CTG or a
        one-base variant); used for the purity filter.
      ttg_fraction: TTG tokens as a fraction of unit-like tokens.
      lower_bound: Repeat units seen in the periodic run (one-anchor reads).
      side: "left" or "right", the anchor that was found (one-anchor reads).
      blocks: Interruption blocks (spanning reads).
      edits: Total anchor edit distance (spanning reads).
    """

    status: str
    reason: str | None = None
    gap: int | None = None
    units: float | None = None
    remainder: int | None = None
    purity: float | None = None
    periodic: float | None = None
    ttg_fraction: float | None = None
    lower_bound: int | None = None
    side: str | None = None
    blocks: tuple[motif.Block, ...] = ()
    edits: int | None = None


def locate_anchor(
    read: str, anchor: str, budget: int
) -> tuple[Hit | None, bool]:
    """Finds the best occurrence of an anchor and tests that it is unique.

    The best hit within ``budget`` edits is found (edlib HW on the whole read,
    O(k n)). It is then masked and the search is repeated; any second hit
    within ``max(budget, best + 3)`` edits makes the anchor ambiguous.

    Args:
      read: Read sequence, uppercase.
      anchor: Anchor sequence, uppercase.
      budget: Maximum edit distance.

    Returns:
      The best hit (None if absent) and True if a second hit exists.
    """
    res = edlib.align(anchor, read, mode="HW", task="locations", k=budget)
    if res["editDistance"] == -1 or not res["locations"]:
        return None, False
    start, end = res["locations"][0]  # edlib end is inclusive
    hit = Hit(start, end + 1, res["editDistance"])
    masked = read[: hit.start] + "N" * (hit.end - hit.start) + read[hit.end :]
    again = edlib.align(
        anchor,
        masked,
        mode="HW",
        task="distance",
        k=max(budget, hit.edits + 3),
    )
    return hit, again["editDistance"] != -1


RESCUE_MAX_EDITS = 25
RESCUE_BEFORE = 50
RESCUE_AFTER = 250


def _rescue(
    seq: str, left: Hit | None, right: Hit | None, run_end: int
) -> tuple[Hit | None, Hit | None]:
    """Looks for the missing anchor as a 100 bp extended anchor.

    Used when one strict anchor is found and its periodic run ends more than
    one window before the read end. The extended anchor shares the inner
    edge of the 40 bp anchor, so the gap is measured from that edge. The
    search is limited to ``[run_end - 50, run_end + 250]`` (read offsets,
    in the direction of the run), allows 25 edits, and applies the same
    uniqueness test as the strict anchors.

    Args:
      seq: Read, plus strand.
      left: Left anchor hit, or None.
      right: Right anchor hit, or None.
      run_end: Read offset where the periodic run ended.

    Returns:
      The (left, right) hits, with the rescued one filled in if found.
    """
    if left is not None and right is None:
        lo = max(0, run_end - RESCUE_BEFORE)
        hit, dup = locate_anchor(
            seq[lo : run_end + RESCUE_AFTER], locus.RIGHT_EXT, RESCUE_MAX_EDITS
        )
        if hit is not None and not dup:
            right = Hit(lo + hit.start, lo + hit.end, hit.edits)
    elif right is not None and left is None:
        lo = max(0, run_end - RESCUE_AFTER)
        hit, dup = locate_anchor(
            seq[lo : run_end + RESCUE_BEFORE], locus.LEFT_EXT, RESCUE_MAX_EDITS
        )
        if hit is not None and not dup:
            left = Hit(lo + hit.start, lo + hit.end, hit.edits)
    return left, right


def measure_read(
    sequence: str,
    max_edits: int = 6,
    min_purity: float = 0.5,
    purity_min_units: float = 20.0,
) -> Measure:
    """Measures the repeat of one read from its own sequence.

    Plus-strand anchors are searched in the sequence as stored (BAM SEQ is
    reference-forward). The reverse complement is searched only to reject
    inverted chimeras and mis-oriented records. If one anchor is missing and
    the periodic run from the other ends before the read does, the missing
    anchor is searched as a 100 bp extended anchor (status "spanning",
    reason "anchor_rescued").

    Args:
      sequence: Primary SEQ, including soft-clipped bases.
      max_edits: Edit budget per strict anchor.
      min_purity: Reads with a lower periodic fraction are rejected ...
      purity_min_units: ... when they have at least this many units.

    Returns:
      A ``Measure``; spanning reads carry gap, units, remainder, purity.
    """
    seq = sequence.upper()
    left, left_dup = locate_anchor(seq, locus.LEFT_ANCHOR, max_edits)
    right, right_dup = locate_anchor(seq, locus.RIGHT_ANCHOR, max_edits)
    if left_dup or right_dup:
        return Measure("rejected", "duplicate_anchor")

    if left is None or right is None:
        rc = locus.reverse_complement(seq)
        rc_left, _ = locate_anchor(rc, locus.LEFT_ANCHOR, max_edits)
        rc_right, _ = locate_anchor(rc, locus.RIGHT_ANCHOR, max_edits)
        if (left is None and rc_left is not None) or (
            right is None and rc_right is not None
        ):
            code = (
                "inverted_chimera"
                if left is not None or right is not None
                else "reverse_complement"
            )
            return Measure("rejected", code)

    if left is None and right is None:
        return Measure("no_anchor", "no_anchor")

    rescued = False
    if left is None or right is None:
        if left is not None:
            origin, direction, side = left.end + locus.MARGIN_LEFT, 1, "left"
            reason = "no_right_anchor"
        else:
            origin, direction, side = (
                right.start - locus.MARGIN_RIGHT,
                -1,
                "right",
            )
            reason = "no_left_anchor"
        covered, to_end = motif.periodic_run(seq, origin, direction)
        if to_end:
            if covered >= 3:
                return Measure(
                    "one_anchor", lower_bound=covered // 3, side=side
                )
            return Measure("partial", reason, side=side)
        left, right = _rescue(seq, left, right, origin + direction * covered)
        if left is None or right is None:
            return Measure("partial", reason, side=side)
        rescued = True

    assert left is not None and right is not None
    gap = right.start - left.end
    if right.start < left.end:
        return Measure("rejected", "order", gap=gap)
    if gap < locus.OFFSET:
        return Measure("rejected", "gap_lt_offset", gap=gap)
    tract = seq[left.end + locus.MARGIN_LEFT : right.start - locus.MARGIN_RIGHT]
    analysis = motif.analyze_tract(tract)
    shifted = gap - locus.OFFSET
    units = shifted / 3
    remainder = shifted % 3
    edits = left.edits + right.edits
    common = {
        "gap": gap,
        "units": units,
        "remainder": remainder,
        "purity": analysis.purity,
        "periodic": analysis.periodic,
        "ttg_fraction": analysis.ttg_fraction,
        "edits": edits,
    }
    if analysis.periodic < min_purity and units >= purity_min_units:
        return Measure("rejected", "low_purity", **common)
    return Measure(
        "spanning",
        reason="anchor_rescued" if rescued else None,
        blocks=analysis.blocks,
        **common,
    )
