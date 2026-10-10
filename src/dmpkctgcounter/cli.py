"""Command line interface: one small TSV plus a summary line on stderr."""

from __future__ import annotations

import argparse
import importlib.metadata
import logging
import os
import platform
import sys
from collections.abc import Sequence

import numpy
import pysam

import dmpkctgcounter
from dmpkctgcounter import locus, pipeline

REFERENCE_ENV = "DMPKCTGCOUNTER_REFERENCE"
COLUMNS = (
    "sample",
    "allele",
    "ctg_call",
    "ctg_median",
    "ctg_min",
    "ctg_max",
    "n_reads",
    "flags",
    "mad",
    "p10",
    "p90",
    "median_fwd",
    "median_rev",
    "n_one_anchor",
    "max_lower_bound",
    "frac_rem0",
    "motif_mode_caa",
)
READ_COLUMNS = (
    "read",
    "strand",
    "mapq",
    "hp",
    "ps",
    "status",
    "reason",
    "gap",
    "units",
    "remainder",
    "purity",
    "periodic",
    "ttg_fraction",
    "lower_bound",
    "side",
    "allele",
)


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dmpkCTGcounter",
        description=(
            "Estimate the DMPK 3' UTR (CTG)n repeat units of both alleles "
            "from one ONT BAM or CRAM (GRCh38; whole-genome or adaptive "
            "sampling). Units = (distance between the two anchors - 7) / 3."
        ),
    )
    p.add_argument("alignment", help="BAM or CRAM with index (GRCh38, ONT)")
    p.add_argument(
        "--reference",
        help="indexed GRCh38 FASTA: verifies the anchor constants and "
        f"decodes CRAM (default: the FASTA named by ${REFERENCE_ENV}, "
        "if set)",
    )
    p.add_argument("--out", help="output TSV (default: stdout)")
    p.add_argument("--sample", help="sample name (default: @RG SM or stem)")
    p.add_argument(
        "--basecall-model",
        help="basecall model name (default: from @RG DS); warns on hac, "
        "fast, R9.4.1, or unknown",
    )
    d = pipeline.Params()
    p.add_argument(
        "--min-mapq",
        type=int,
        default=d.min_mapq,
        help="skip primary records below this MAPQ (default %(default)s)",
    )
    p.add_argument(
        "--max-edits",
        type=int,
        default=d.max_edits,
        help="edit budget per anchor (default %(default)s)",
    )
    p.add_argument(
        "--min-purity",
        type=float,
        default=d.min_purity,
        help="reject reads whose unit-like token fraction is below this "
        "if they have at least "
        "--purity-min-units units (default %(default)s)",
    )
    p.add_argument(
        "--purity-min-units",
        type=float,
        default=d.purity_min_units,
        help="see --min-purity (default %(default)s)",
    )
    p.add_argument(
        "--per-read",
        metavar="TSV",
        help="also write one row per read to this file",
    )
    p.add_argument(
        "--size-class",
        action="store_true",
        help="add a size class column (research use only)",
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {dmpkctgcounter.__version__}",
    )
    return p


def _num(value: float | None, digits: int = 1) -> str:
    return "NA" if value is None else f"{value:.{digits}f}"


def size_class(
    ctg_call: int | None, n_reads: int = 0, median: float | None = None
) -> str:
    """ACMG 2024 style size class of a repeat count (research use only).

    NA when there is no call, or when the allele rests on one read below 50
    units (a single short read is not a size estimate).
    """
    if ctg_call is None:
        return "NA"
    if n_reads == 1 and (median is None or median < 50):
        return "NA"
    if ctg_call < 5:
        return "<5"
    if ctg_call <= 34:
        return "normal(5-34)"
    if ctg_call <= 49:
        return "variable(35-49)"
    return "fully_variable(>=50)"


def format_tsv(
    result: pipeline.RunResult,
    args: argparse.Namespace,
    argv: Sequence[str],
) -> str:
    """Renders the TSV: `#` header lines, a column line, one row per allele."""
    lines = [
        f"# dmpkCTGcounter {dmpkctgcounter.__version__}",
        f"# command: dmpkCTGcounter {' '.join(argv)}",
        f"# versions: python {platform.python_version()}; "
        f"pysam {pysam.__version__}; "
        f"edlib {importlib.metadata.version('edlib')}; "
        f"numpy {numpy.__version__}",
        f"# input: {args.alignment} ({result.input_format}); "
        f"reference: {result.reference_status}"
        + (f" [from {REFERENCE_ENV}]" if args.reference_auto else ""),
        f"# build: GRCh38 (chr19 length {locus.GRCH38_CHR19_LENGTH} "
        f"asserted; alignment contig {result.alignment_contig})",
        f"# tract: {locus.TRACT_HGVS}; reference units {locus.REF_UNITS}",
        f"# left anchor chr19:{locus.LEFT_START}-{locus.LEFT_END} "
        f"{locus.LEFT_ANCHOR}; right anchor "
        f"chr19:{locus.RIGHT_START}-{locus.RIGHT_END} {locus.RIGHT_ANCHOR}",
        f"# units = (gap - {locus.OFFSET}) / 3; reference gap "
        f"{locus.REF_GAP} bp = {(locus.REF_GAP - locus.OFFSET) / 3:.1f} "
        "units; ctg_call = floor(median + 0.5); no bias correction",
        f"# basecall model: {result.basecall_model or 'unknown'}",
        f"# params: min_mapq={args.min_mapq}; max_edits={args.max_edits}; "
        f"min_purity={args.min_purity}; "
        f"purity_min_units={args.purity_min_units}",
        f"# method: {result.method}; deterministic, no random seed",
        "# counts: "
        + "; ".join(f"{k}={v}" for k, v in sorted(result.counts.items()))
        + (
            "; reasons: "
            + ", ".join(f"{k}={v}" for k, v in sorted(result.reasons.items()))
            if result.reasons
            else ""
        ),
    ]
    lines += [f"# warning: {w}" for w in result.warnings]
    columns = list(COLUMNS)
    if args.size_class:
        columns.append("size_class_research_use")
    lines.append("\t".join(columns))
    for row in result.rows:
        fields = [
            result.sample,
            str(row.allele),
            "NA" if row.ctg_call is None else str(row.ctg_call),
            _num(row.median),
            _num(row.minimum),
            _num(row.maximum),
            str(row.n_reads),
            ";".join(row.flags) if row.flags else ".",
            _num(row.mad, 2),
            _num(row.p10),
            _num(row.p90),
            _num(row.median_fwd),
            _num(row.median_rev),
            str(result.n_one_anchor),
            "NA"
            if result.max_lower_bound is None
            else str(result.max_lower_bound),
            _num(row.frac_rem0, 2),
            "NA" if row.motif_mode_caa is None else str(row.motif_mode_caa),
        ]
        if args.size_class:
            fields.append(size_class(row.ctg_call, row.n_reads, row.median))
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def format_per_read(result: pipeline.RunResult) -> str:
    """Renders the optional per-read table."""
    lines = ["\t".join(READ_COLUMNS)]
    for r in result.reads:
        m = r.measure
        lines.append(
            "\t".join(
                [
                    r.name,
                    r.strand,
                    str(r.mapq),
                    "NA" if r.hp is None else str(r.hp),
                    "NA" if r.ps is None else str(r.ps),
                    m.status,
                    m.reason or ".",
                    "NA" if m.gap is None else str(m.gap),
                    _num(m.units, 2),
                    "NA" if m.remainder is None else str(m.remainder),
                    _num(m.purity, 2),
                    _num(m.periodic, 2),
                    _num(m.ttg_fraction, 2),
                    "NA" if m.lower_bound is None else str(m.lower_bound),
                    m.side or ".",
                    "NA" if r.allele is None else str(r.allele),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def summary_line(result: pipeline.RunResult) -> str:
    """One human-readable line for stderr."""
    parts = []
    for row in result.rows:
        if row.median is None:
            parts.append(f"allele{row.allele} NA [{';'.join(row.flags)}]")
        else:
            parts.append(
                f"allele{row.allele} {row.ctg_call} CTG (median "
                f"{row.median:.1f}, range {row.minimum:.1f}-"
                f"{row.maximum:.1f}, n={row.n_reads})"
                + (f" [{';'.join(row.flags)}]" if row.flags else "")
            )
    c = result.counts
    text = (
        f"{result.sample}: "
        + "; ".join(parts)
        + f"; spanning reads {c.get('spanning', 0)} of {c.get('reads', 0)}"
        + f"; one-anchor reads {result.n_one_anchor}"
    )
    if result.max_lower_bound:
        text += f" (longest run >= {result.max_lower_bound} units)"
    return text


def main(argv: Sequence[str] | None = None) -> int:
    """Runs the command line tool; returns the process exit code."""
    argv = list(sys.argv[1:] if argv is None else argv)
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    reference = args.reference
    args.reference_auto = False
    env_reference = os.environ.get(REFERENCE_ENV)
    if reference is None and env_reference:
        if not os.path.exists(env_reference + ".fai"):
            print(
                f"dmpkCTGcounter: error: {REFERENCE_ENV} is set to "
                f"{env_reference!r}, but its .fai index was not found",
                file=sys.stderr,
            )
            return 2
        reference = env_reference
        args.reference_auto = True
    args.reference = reference
    params = pipeline.Params(
        min_mapq=args.min_mapq,
        max_edits=args.max_edits,
        min_purity=args.min_purity,
        purity_min_units=args.purity_min_units,
    )
    try:
        result = pipeline.count_ctg(
            args.alignment,
            reference,
            params,
            args.sample,
            args.basecall_model,
        )
    except locus.InputError as err:
        print(f"dmpkCTGcounter: error: {err}", file=sys.stderr)
        return 2
    text = format_tsv(result, args, argv)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    if args.per_read:
        with open(args.per_read, "w", encoding="utf-8") as handle:
            handle.write(format_per_read(result))
    print(summary_line(result), file=sys.stderr)
    for warning in result.warnings:
        print(f"dmpkCTGcounter: warning: {warning}", file=sys.stderr)
    return 0
