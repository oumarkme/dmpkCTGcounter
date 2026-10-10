"""DMPK (CTG)n locus constants for GRCh38 and build checks.

Coordinate conventions
----------------------
* All genomic positions in this module and in the output are 1-based and
  closed (SAM and VCF convention), on ``chr19`` of GRCh38.
* Offsets into read strings elsewhere in the package are 0-based half-open;
  ``verify_reference`` is the one place where a FASTA is sliced, using
  ``pysam.FastaFile.fetch`` (0-based half-open) with an explicit conversion.
* DMPK lies on the minus strand. The reference plus strand carries the tract
  as (CAG)n, the gene strand as (CTG)n. Repeat units are reported as CTG.

Reference structure around the tract (plus strand, GRCh38)::

    left anchor 45,770,160-45,770,199 (40 bp; extended 45,770,100-45,770,199)
    margin      45,770,200-45,770,204  CCCCC            (5 bp)
    tract       45,770,205-45,770,264  (CAG)20          (60 bp)
    margin      45,770,265-45,770,266  CA               (2 bp)
    right anchor 45,770,267-45,770,306 (40 bp; extended 45,770,267-45,770,366)

The gap between the inner anchor edges is 67 bp = 7 bp of margin + 3n bp of
repeat, so the repeat count of a read is ``(gap - 7) / 3``.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping

import pysam

LOG = logging.getLogger("dmpkctgcounter")

GRCH38_CHR19_LENGTH = 58_617_616
GRCH38_CHR1_LENGTH = 248_956_422
# Explicit contig map (no stripping of "chr"): UCSC and Ensembl/NCBI names.
CHR19_NAMES = ("chr19", "19")
CHR1_NAMES = ("chr1", "1")

LEFT_ANCHOR = "AGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGAT"
LEFT_START = 45_770_160
LEFT_END = 45_770_199
RIGHT_ANCHOR = "TTCCCGGCTACAAGGACCCTTCGAGCCCCGTTCGCCGGCC"
RIGHT_START = 45_770_267
RIGHT_END = 45_770_306
# Extended anchors (100 bp, same inner edges) used only to rescue reads whose
# 40 bp anchor is degraded; checked against the reference like the others.
LEFT_EXT_START = 45_770_100
LEFT_EXT = (
    "ATGGAACACGGACGGCCCGGCTTGCTGCCTTCCCAGGCCTGCAGTTTGCC"
    "CATCCACGTCAGGGCCTCAGCCTGGCCGAAAGAAAGAAATGGTCTGTGAT"
)
RIGHT_EXT_END = 45_770_366
RIGHT_EXT = (
    "TTCCCGGCTACAAGGACCCTTCGAGCCCCGTTCGCCGGCCGCGGACCCGG"
    "CCCCTCCCTCCCCGGCCGCTAGGGGGCGGGCCCGGATCACAGGACTGGAG"
)
TRACT_START = 45_770_205
TRACT_END = 45_770_264
REF_UNITS = 20
MARGIN_LEFT = 5  # bp between left anchor and tract (CCCCC)
MARGIN_RIGHT = 2  # bp between tract and right anchor (CA)
OFFSET = MARGIN_LEFT + MARGIN_RIGHT
REF_GAP = RIGHT_START - LEFT_END - 1  # 67
PLUS_MOTIF = "CAG"
GENE_MOTIF = "CTG"
TRACT_HGVS = "NC_000019.10:g.45770205_45770264 = NM_004409.5:c.*224_*283CTG[20]"
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class InputError(Exception):
    """A user-correctable problem with the inputs (reported, exit code 2)."""


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of an uppercase ACGTN string."""
    return seq.translate(_COMPLEMENT)[::-1]


def resolve_contig(
    names_to_lengths: Mapping[str, int], aliases: tuple[str, ...], what: str
) -> str | None:
    """Returns the single contig name among ``aliases``, or None if absent.

    Args:
      names_to_lengths: Contig name to length, from a header or FASTA index.
      aliases: Accepted names for the contig (explicit map).
      what: Label used in error messages, e.g. "chr19".

    Raises:
      InputError: If more than one alias is present.
    """
    found = [name for name in aliases if name in names_to_lengths]
    if len(found) > 1:
        raise InputError(f"both {found} present; cannot resolve {what}")
    return found[0] if found else None


def check_build(names_to_lengths: Mapping[str, int], source: str) -> str:
    """Asserts GRCh38 from contig lengths and returns the chr19 name.

    Args:
      names_to_lengths: Contig name to length (BAM/CRAM @SQ or FASTA index).
      source: Description of where the lengths come from, for messages.

    Raises:
      InputError: If chr19 is missing or its length (or chr1's, when present)
        is not the GRCh38 value.
    """
    contig = resolve_contig(names_to_lengths, CHR19_NAMES, "chr19")
    if contig is None:
        raise InputError(
            f"{source}: no contig named any of {CHR19_NAMES}; "
            "GRCh38 chr19 is required"
        )
    length = names_to_lengths[contig]
    if length != GRCH38_CHR19_LENGTH:
        raise InputError(
            f"{source}: {contig} has length {length}, expected "
            f"{GRCH38_CHR19_LENGTH} (GRCh38); other builds are not supported"
        )
    chr1 = resolve_contig(names_to_lengths, CHR1_NAMES, "chr1")
    if chr1 is not None and names_to_lengths[chr1] != GRCH38_CHR1_LENGTH:
        raise InputError(
            f"{source}: {chr1} has length {names_to_lengths[chr1]}, expected "
            f"{GRCH38_CHR1_LENGTH} (GRCh38)"
        )
    return contig


def verify_reference(reference: str) -> str:
    """Checks that a FASTA carries the anchor constants and the tract.

    Args:
      reference: Path to an indexed FASTA.

    Returns:
      The chr19 contig name in the FASTA.

    Raises:
      InputError: If the build is not GRCh38 or any base differs from the
        constants (anchors, margins, or the (CAG)20 tract).
    """
    with pysam.FastaFile(reference) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        contig = check_build(lengths, f"reference {reference}")
        # 1-based closed LEFT_EXT_START..RIGHT_EXT_END -> 0-based half-open.
        seq = fasta.fetch(contig, LEFT_EXT_START - 1, RIGHT_EXT_END).upper()
    base = LEFT_EXT_START
    expected = [
        ("left anchor", LEFT_START, LEFT_ANCHOR),
        ("left margin", TRACT_START - MARGIN_LEFT, "CCCCC"),
        ("tract", TRACT_START, PLUS_MOTIF * REF_UNITS),
        ("right margin", TRACT_END + 1, "CA"),
        ("right anchor", RIGHT_START, RIGHT_ANCHOR),
        ("left extended anchor", LEFT_EXT_START, LEFT_EXT),
        ("right extended anchor", RIGHT_START, RIGHT_EXT),
    ]
    for label, start, bases in expected:
        found = seq[start - base : start - base + len(bases)]
        if found != bases:
            first = next(
                i
                for i, (a, b) in enumerate(zip(found, bases, strict=False))
                if a != b
            )
            raise InputError(
                f"reference {reference}: {label} differs from the GRCh38 "
                f"sequence at {contig}:{start + first} (found "
                f"{found[first]!r}, expected {bases[first]!r}); this is not "
                "the GRCh38 reference the anchors were defined on"
            )
    return contig
