"""BAM/CRAM validation, read selection, and the end-to-end run.

Read selection (one region fetch, no whole-file scan):

* One fetch of chr19:45,760,000-45,780,000 (1-based closed; converted once
  to the 0-based half-open interval ``(45759999, 45780000)``).
* Records with flags 0x4, 0x100, 0x200, or 0x400 are dropped, and so are
  duplex parent reads (``dx`` tag equal to -1).
* Reads are keyed by name. The primary record's full SEQ (soft-clipped bases
  included) is measured. A supplementary record found in the region is
  resolved to its primary through the SA tag; if the two lie on opposite
  strands, the primary SEQ is reverse-complemented so that it is
  reference-forward at chr19. A primary with hard clips has no full SEQ and
  is skipped and counted.
* HP tags are used only within the dominant PS phase block among the
  spanning reads; reads of another block are treated as untagged.
"""

from __future__ import annotations

import collections
import dataclasses
import hashlib
import logging
import os
import re

import pysam

from dmpkctgcounter import calling, locus, matching

LOG = logging.getLogger("dmpkctgcounter")

FETCH_START = 45_760_000  # 1-based closed
FETCH_END = 45_780_000
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400
_SUPPLEMENTARY = 0x800
_MAX_SA_LOOKUPS = 20
_HARD_CLIP = 5  # BAM CIGAR operation H
_MODEL_FULL = re.compile(
    r"dna_r\d+\.\d+\.\d+_e\d+\.\d+_\d+bps_(?:fast|hac|sup)@v\d+\.\d+\.\d+"
)
_MODEL_SHORT = re.compile(r"^(fast|hac|sup)(@v\d+\.\d+\.\d+)?$")
_MODBASE = re.compile(r"_(?:5mC|5hmC|6mA|4mC|5mCG|m\d)")


@dataclasses.dataclass(frozen=True)
class Params:
    """Run parameters (all are logged in the TSV header)."""

    min_mapq: int = 0
    max_edits: int = 6
    min_purity: float = 0.5
    purity_min_units: float = 20.0


@dataclasses.dataclass
class ReadRow:
    """Per-read result (for ``--per-read``)."""

    name: str
    strand: str
    mapq: int
    hp: int | None
    ps: int | None
    measure: matching.Measure
    allele: int | None = None


@dataclasses.dataclass
class RunResult:
    """Everything the writers need."""

    sample: str
    input_format: str
    alignment_contig: str
    reference_status: str
    basecall_model: str | None
    rows: list[calling.AlleleCall]
    method: str
    counts: dict[str, int]
    reasons: collections.Counter
    reads: list[ReadRow]
    n_one_anchor: int
    max_lower_bound: int | None
    warnings: list[str]


def sniff_format(path: str) -> str:
    """Returns "BAM" or "CRAM" from the first bytes of the file.

    Raises:
      locus.InputError: If the file is unreadable or neither BAM nor CRAM.
    """
    try:
        with open(path, "rb") as handle:
            magic = handle.read(4)
    except OSError as err:
        raise locus.InputError(f"cannot read {path}: {err}") from err
    if magic == b"CRAM":
        return "CRAM"
    if magic[:2] == b"\x1f\x8b":
        return "BAM"
    raise locus.InputError(f"{path} is not a BAM or CRAM file")


def _sample_name(header: pysam.AlignmentHeader, path: str) -> str:
    """Returns the single @RG SM, or the file name without extension."""
    names = {rg.get("SM") for rg in header.to_dict().get("RG", [])}
    names.discard(None)
    if len(names) == 1:
        return names.pop()
    base = os.path.basename(path)
    for suffix in (".bam", ".cram"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def _basecall_model(header: pysam.AlignmentHeader) -> str | None:
    """Collects the simplex basecall model names of the header.

    Sources, all of them used: every dorado @PG command line (a
    ``dna_r..._sup@vX.Y.Z`` token, skipping modified-base model names, with the
    dorado version appended; a bare ``sup`` or ``hac`` token is kept as such),
    and ``basecall_model=`` in every @RG DS. If nothing is found the result is
    "unknown (MinKNOW <version>)" from the MinKNOW @PG entry.

    Args:
      header: Alignment header.

    Returns:
      Comma-joined model names, or None if the header has no relevant field.
    """
    meta = header.to_dict()
    models: list[str] = []
    minknow = None
    for pg in meta.get("PG", []):
        name = str(pg.get("PN", "")).lower()
        if name == "minknow" and minknow is None:
            minknow = pg.get("VN")
        if name != "dorado":
            continue
        version = pg.get("VN", "unknown")
        for token in str(pg.get("CL", "")).split():
            token = os.path.basename(token.strip("'\""))
            if _MODBASE.search(token):
                continue
            full = _MODEL_FULL.search(token)
            if full:
                models.append(f"{full.group(0)} (dorado {version})")
                break
            if _MODEL_SHORT.match(token):
                models.append(f"{token} (dorado {version})")
                break
    for rg in meta.get("RG", []):
        found = re.search(r"basecall_model=(\S+)", rg.get("DS", ""))
        if found:
            models.append(found.group(1))
    models = list(dict.fromkeys(models))
    if models:
        return ",".join(models)
    if minknow is not None:
        return f"unknown (MinKNOW {minknow})"
    return None


def model_warning(model: str | None) -> str | None:
    """Returns a warning if any basecall model may bias repeat counts."""
    if not model or model.startswith("unknown"):
        return (
            f"basecall model {model or 'unknown'}; use --basecall-model. "
            "Counts from fast or hac models, or R9.4.1 chemistry, may be "
            "biased"
        )
    notes = []
    for part in model.split(","):
        low = part.lower()
        if "dna_r" not in low:
            notes.append(
                f"basecall model {part} names no chemistry (chemistry "
                "unknown); R9.4.1 or non-sup counts may be biased"
            )
        elif "sup" not in low or "r9" in low or "9.4.1" in low:
            notes.append(
                f"basecall model {part} is not an R10 sup model; repeat "
                "counts may be biased (indel errors in the repeat)"
            )
        else:
            version = re.search(r"@v(\d+)\.", low)
            if version and int(version.group(1)) < 5:
                notes.append(
                    f"basecall model {part} is a sup model earlier than v5; "
                    "expanded repeats may be called CAA-rich (CTG to TTG) "
                    "and single-read counts may be biased"
                )
    return "; ".join(notes) if notes else None


def _reference_md5(reference: str, contig: str) -> str:
    """MD5 of the uppercase contig sequence, as in the SAM @SQ M5 field."""
    with pysam.FastaFile(reference) as fasta:
        seq = fasta.fetch(contig).upper()
    return hashlib.md5(seq.encode("ascii")).hexdigest()


def _primary_from_sa(
    alignment: pysam.AlignmentFile, record: pysam.AlignedSegment
) -> pysam.AlignedSegment | None:
    """Follows the SA tag of a supplementary record to the primary record."""
    if not record.has_tag("SA"):
        return None
    entries = [e for e in str(record.get_tag("SA")).split(";") if e]
    for entry in entries[:_MAX_SA_LOOKUPS]:
        fields = entry.split(",")
        if len(fields) < 4 or fields[0] not in alignment.references:
            continue
        pos0 = int(fields[1]) - 1  # SA POS is 1-based
        for cand in alignment.fetch(fields[0], pos0, pos0 + 1):
            if cand.query_name == record.query_name and not (
                cand.flag & (_SKIP_FLAGS | _SUPPLEMENTARY)
            ):
                return cand
    return None


def _tag(rec: pysam.AlignedSegment, name: str) -> int | None:
    """Returns an integer tag value or None."""
    if not rec.has_tag(name):
        return None
    return int(rec.get_tag(name))


@dataclasses.dataclass
class _Candidate:
    """A read to measure: the primary record and how to orient its SEQ."""

    primary: pysam.AlignedSegment
    flip: bool
    reverse: bool


def count_ctg(
    path: str,
    reference: str | None = None,
    params: Params | None = None,
    sample: str | None = None,
    basecall_model: str | None = None,
) -> RunResult:
    """Estimates the CTG repeat units of the two DMPK alleles.

    Args:
      path: BAM or CRAM with an index (.bai, .csi, or .crai); GRCh38; ONT.
      reference: Indexed GRCh38 FASTA. Needed to decode CRAM unless the CRAM
        embeds its reference; if given, the anchor constants are verified
        against it. Without it the built-in constants are used.
      params: Run parameters; defaults if None.
      sample: Optional sample name; default is the @RG SM or the file stem.
      basecall_model: Optional model name; default is read from the header.

    Returns:
      A ``RunResult``.

    Raises:
      locus.InputError: Wrong build, missing index or reference, bad file,
        or a record that cannot be decoded.
    """
    params = params or Params()
    fmt = sniff_format(path)
    warnings: list[str] = []
    if reference is not None:
        if not os.path.exists(reference) or not os.path.exists(
            reference + ".fai"
        ):
            raise locus.InputError(
                f"reference {reference} not found or has no .fai index "
                "(run `samtools faidx`)"
            )
        ref_contig = locus.verify_reference(reference)
        ref_status = f"{reference} (anchor constants verified)"
    else:
        ref_contig = None
        ref_status = "built-in constants (no reference checked)"

    try:
        alignment = pysam.AlignmentFile(
            path, "rb", reference_filename=reference
        )
    except (OSError, ValueError) as err:
        raise locus.InputError(f"cannot open {path}: {err}") from err
    candidates: dict[str, _Candidate] = {}
    counts: collections.Counter = collections.Counter()
    with alignment:
        if not alignment.has_index():
            raise locus.InputError(
                f"no index for {path}; run `samtools index {path}`"
            )
        lengths = dict(
            zip(alignment.references, alignment.lengths, strict=True)
        )
        contig = locus.check_build(lengths, f"alignment {path}")
        sq_entries = {
            e["SN"]: e for e in alignment.header.to_dict().get("SQ", [])
        }
        sample_name = sample or _sample_name(alignment.header, path)
        model = basecall_model or _basecall_model(alignment.header)
        note = model_warning(model)
        if note:
            warnings.append(note)
        m5 = sq_entries[contig].get("M5")
        if (
            m5
            and reference is not None
            and m5 != _reference_md5(reference, ref_contig)
        ):
            warnings.append(
                f"@SQ M5 of {contig} differs from the chr19 sequence of "
                "--reference; the alignment was made against another FASTA"
            )

        primaries: dict[str, pysam.AlignedSegment] = {}
        supplementary: dict[str, pysam.AlignedSegment] = {}
        try:
            for rec in alignment.fetch(contig, FETCH_START - 1, FETCH_END):
                counts["records_fetched"] += 1
                if rec.flag & _SKIP_FLAGS:
                    counts["flag_skipped_records"] += 1
                elif _tag(rec, "dx") == -1:
                    counts["duplex_parent"] += 1
                elif rec.flag & _SUPPLEMENTARY:
                    supplementary.setdefault(rec.query_name, rec)
                else:
                    primaries.setdefault(rec.query_name, rec)
            for name, rec in primaries.items():
                candidates[name] = _Candidate(rec, False, rec.is_reverse)
            for name, srec in supplementary.items():
                if name in candidates:
                    continue
                prim = _primary_from_sa(alignment, srec)
                if prim is None:
                    counts["primary_not_found"] += 1
                elif _tag(prim, "dx") == -1:
                    counts["duplex_parent"] += 1
                else:
                    candidates[name] = _Candidate(
                        prim,
                        prim.is_reverse != srec.is_reverse,
                        srec.is_reverse,
                    )
        except (OSError, ValueError) as err:
            if fmt == "CRAM" and reference is None:
                raise locus.InputError(
                    f"cannot decode CRAM {path} without its reference "
                    f"({err}); give --reference"
                ) from err
            if fmt == "CRAM":
                raise locus.InputError(
                    f"cannot decode CRAM {path} with reference {reference} "
                    f"({err}); the FASTA may differ from the one used to "
                    "encode it (MD5 mismatch)"
                ) from err
            raise locus.InputError(
                f"cannot decode records of {path} ({err})"
            ) from err

    reasons: collections.Counter = collections.Counter()
    reads: list[ReadRow] = []
    for name, cand in candidates.items():
        rec = cand.primary
        cigar = rec.cigartuples or []
        seq = rec.query_sequence
        if any(op == _HARD_CLIP for op, _ in cigar) or not seq:
            counts["hard_clipped_primary"] += 1
            continue
        if rec.mapping_quality < params.min_mapq:
            counts["below_min_mapq"] += 1
            continue
        if rec.mapping_quality < 20:
            counts["mapq_lt_20"] += 1
        if cand.flip:
            seq = locus.reverse_complement(seq.upper())
        measure = matching.measure_read(
            seq,
            params.max_edits,
            params.min_purity,
            params.purity_min_units,
        )
        reads.append(
            ReadRow(
                name,
                "-" if cand.reverse else "+",
                rec.mapping_quality,
                _tag(rec, "HP"),
                _tag(rec, "PS"),
                measure,
            )
        )
        counts[measure.status] += 1
        if measure.reason:
            reasons[measure.reason] += 1
    counts["reads"] = len(reads)

    spanning = [r for r in reads if r.measure.status == "spanning"]
    phased = [r for r in spanning if r.hp in (1, 2)]
    dominant = None
    if phased:
        tally = collections.Counter(r.ps for r in phased)
        dominant = sorted(tally, key=lambda k: (-tally[k], str(k)))[0]
    points = []
    for r in spanning:
        m = r.measure
        hp = r.hp if r.hp in (1, 2) and r.ps == dominant else None
        if r.hp in (1, 2) and hp is None:
            counts["hp_other_phase_block"] += 1
        points.append(
            calling.ReadPoint(
                r.name,
                m.units,
                hp,
                r.strand == "-",
                m.periodic or 0.0,
                m.remainder or 0,
                m.blocks,
                strict_purity=m.purity or 0.0,
                ttg_fraction=m.ttg_fraction or 0.0,
                rescued=m.reason == "anchor_rescued",
            )
        )
    bounds = [
        r.measure.lower_bound
        for r in reads
        if r.measure.status == "one_anchor" and r.measure.lower_bound
    ]
    max_lb = max(bounds) if bounds else None
    rows, method, _ = calling.call_alleles(points, max_lb)
    if not reads:
        marker = (
            "NO_READS_IN_REGION"
            if counts["records_fetched"] == 0
            else "NO_USABLE_READS"
        )
        for row in rows:
            row.flags = [marker]
    by_name = {p.name: p.allele for p in points}
    for r in reads:
        r.allele = by_name.get(r.name)
    return RunResult(
        sample=sample_name,
        input_format=fmt,
        alignment_contig=contig,
        reference_status=ref_status,
        basecall_model=model,
        rows=rows,
        method=method,
        counts=dict(counts),
        reasons=reasons,
        reads=reads,
        n_one_anchor=len(bounds),
        max_lower_bound=max_lb,
        warnings=warnings,
    )
