"""Allele calling from per-read repeat counts of spanning reads.

Input values are real repeat counts ``u = (gap - 7) / 3``. Calls use the
median of ``u`` per allele; no trimming and no bias correction are applied.
Grouping and acceptance work in log space (``log(u + 1)``) because counting
noise grows with the repeat length. The procedure is deterministic
(exhaustive split search, no random numbers). Thresholds are expressed in
reads, units, and log ratios, so the same rules apply to whole-genome and
adaptive-sampling data.
"""

from __future__ import annotations

import collections
import dataclasses
import math
import statistics
from collections.abc import Callable, Sequence

import numpy

from dmpkctgcounter import motif

MAD_SCALE = 1.4826
MIN_HP_FRACTION = 0.6
MIN_GROUP_READS = 2
# Two alleles must differ by at least this many units (at the lower median,
# expressed as a log ratio) and by this many pooled robust SDs in log space.
MIN_DELTA_UNITS = 2.0
DELTA_MAD_MULTIPLE = 3.0
SINGLE_READ_MIN_UNITS = 50.0
SINGLE_READ_MIN_PURITY = 0.8
HP_DISCORDANT_FRACTION = 0.2
HP_OVERRIDE_DELTA = 10.0
HP_OVERRIDE_RATIO = 2.0
HP_CLOSE_MIN_READS = 3
LOW_DEPTH_READS = 5
LONGER_UNSPANNED_MARGIN = 5
MIXED_RATIO = 2.0
MIXED_RANGE_UNITS = 10.0
EXPANDED_MIN_UNITS = 50.0
EXPANDED_LOW_SUPPORT_READS = 5
STRAND_MIN_READS = 3
STRAND_MIN_LOG = 0.02
STRAND_SE_MULTIPLE = 3.0
P_PERCENTILE_MIN_READS = 10
CANDIDATE_MIN_READS = 3
CANDIDATE_MIN_FRACTION = 0.5
CANDIDATE_OFFSET_TOLERANCE = 3
CANDIDATE_MIN_PURITY = 0.5
CAA_MODE_TTG_FRACTION = 0.5
# Splits of two normal-size groups (median below 50 units) must also look like
# a real heterozygote and not like a noisy single allele cut in two: the
# smaller group must hold at least 25% of the reads, and the medians must
# differ by this many pooled within-group SDs (classical SD, units). The
# threshold falls with the number of reads because the chance of a large
# apparent split in one allele falls with it (seeded simulations).
MIN_ALLELE_FRACTION = 0.25
DPRIME_BY_READS = ((12, 6.0), (20, 5.5), (40, 5.0))  # (n below, threshold)
DPRIME_DEFAULT = 4.5
SD_FLOOR_UNITS = 0.01
# A failed split of two normal-size groups closer than this is "close"
# (and a cut closer than CLOSE_NOISE_UNITS is not reported at all).
CLOSE_MAX_UNITS = 5.0
CLOSE_NOISE_UNITS = 0.3
# An expanded cluster (median 50+) this many times longer than a normal allele
# is a separate allele whatever the read counts.
CLEAR_RATIO = 1.5
# A class inside an already accepted allele is split off only with this many
# reads (and at least CLOSE_MAX_UNITS apart).
RECUT_MIN_READS = 3
# A minor group is a cluster (not a scatter) if its range is at most this.
TIGHT_MIN_RANGE_UNITS = 3.0
TIGHT_RANGE_FRACTION = 0.25
_EPS = 1e-9


@dataclasses.dataclass
class ReadPoint:
    """A spanning read that passed QC.

    Attributes:
      name: Read name.
      units: Repeat count, real.
      hp: HP tag (1 or 2) within the dominant phase block, else None.
      reverse: True if the primary alignment is on the reverse strand.
      purity: Fraction of the tract in unit-like tokens (CTG or a one-base
        variant).
      remainder: ``(gap - 7) mod 3``.
      blocks: Interruption blocks.
      strict_purity: Strict CTG token fraction.
      ttg_fraction: TTG tokens as a fraction of unit-like tokens.
      rescued: True if one anchor came from the extended-anchor rescue.
      allele: Assigned allele (1 or 2), set by ``call_alleles``.
    """

    name: str
    units: float
    hp: int | None
    reverse: bool
    purity: float
    remainder: int
    blocks: tuple[motif.Block, ...] = ()
    strict_purity: float = 1.0
    ttg_fraction: float = 0.0
    rescued: bool = False
    allele: int | None = None


@dataclasses.dataclass
class AlleleCall:
    """One allele row; ``None`` fields are written as NA."""

    allele: int
    ctg_call: int | None = None
    median: float | None = None
    minimum: float | None = None
    maximum: float | None = None
    n_reads: int = 0
    mad: float | None = None
    p10: float | None = None
    p90: float | None = None
    median_fwd: float | None = None
    median_rev: float | None = None
    frac_rem0: float | None = None
    motif_mode_caa: int | None = None
    flags: list[str] = dataclasses.field(default_factory=list)


def _y(units: float) -> float:
    """Log-space value used for grouping."""
    return math.log(max(units, 0.0) + 1.0)


def _mad(values: Sequence[float]) -> float:
    med = statistics.median(values)
    return statistics.median([abs(v - med) for v in values])


def _median_units(group: Sequence[ReadPoint]) -> float:
    return statistics.median(p.units for p in group)


def _accepts_log(low: Sequence[ReadPoint], high: Sequence[ReadPoint]) -> bool:
    """Tests that two groups are distinct in log space.

    The log difference of the medians must be at least the log ratio of a
    2-unit step at the lower median, and at least three pooled robust SDs
    of the within-group log values (MAD x 1.4826). The first term is the
    floor for exact reads; the second keeps a broad single allele whole.
    """
    med_low = _median_units(low)
    med_high = _median_units(high)
    delta = _y(med_high) - _y(med_low)
    pooled = math.sqrt(
        (
            _mad([_y(p.units) for p in low]) ** 2
            + _mad([_y(p.units) for p in high]) ** 2
        )
        / 2
    )
    floor = math.log((med_low + 1 + MIN_DELTA_UNITS) / (med_low + 1))
    return delta >= max(floor, DELTA_MAD_MULTIPLE * MAD_SCALE * pooled) - _EPS


def _dprime_threshold(n_reads: int) -> float:
    """Separation (in pooled SDs) required for a normal-size split."""
    for limit, value in DPRIME_BY_READS:
        if n_reads < limit:
            return value
    return DPRIME_DEFAULT


def _best_cut(
    pts: Sequence[ReadPoint], allowed: Callable[[int], bool]
) -> int | None:
    """Exhaustive 2-means cut of sorted ``log(u + 1)``; ties go lowest.

    Args:
      pts: Reads sorted by ascending units.
      allowed: Predicate on the cut index ``i`` (lower ``[:i]``, upper
        ``[i:]``).

    Returns:
      The cut index with the smallest within-group sum of squares, or None.
    """
    n = len(pts)
    csum = [0.0]
    csq = [0.0]
    for p in pts:
        y = _y(p.units)
        csum.append(csum[-1] + y)
        csq.append(csq[-1] + y * y)

    def sse(lo: int, hi: int) -> float:
        total = csum[hi] - csum[lo]
        return csq[hi] - csq[lo] - total * total / (hi - lo)

    best_i = None
    best_cost = math.inf
    for i in range(1, n):
        if not allowed(i):
            continue
        cost = sse(0, i) + sse(i, n)
        if cost < best_cost - 1e-12:
            best_cost = cost
            best_i = i
    return best_i


def _pooled_sd(low: Sequence[ReadPoint], high: Sequence[ReadPoint]) -> float:
    """Classical pooled within-group SD in units (0 if no group has 2 reads)."""
    ss = 0.0
    dof = 0
    for group in (low, high):
        if len(group) > 1:
            values = [p.units for p in group]
            ss += statistics.variance(values) * (len(values) - 1)
            dof += len(values) - 1
    return math.sqrt(ss / dof) if dof else 0.0


def _is_tight(group: Sequence[ReadPoint]) -> bool:
    """True if a small group looks like one cluster, not a scatter."""
    values = [p.units for p in group]
    median = statistics.median(values)
    return max(values) - min(values) <= max(
        TIGHT_MIN_RANGE_UNITS, TIGHT_RANGE_FRACTION * median
    )


@dataclasses.dataclass
class _Notes:
    """Reads set aside and cuts rejected while segmenting."""

    outliers: list[ReadPoint] = dataclasses.field(default_factory=list)
    named: list[list[ReadPoint]] = dataclasses.field(default_factory=list)
    flags: list[str] = dataclasses.field(default_factory=list)


def _far_from(read: ReadPoint, rest: Sequence[ReadPoint]) -> bool:
    """True if a single read is 10+ units and a factor of 2 from the rest."""
    med = _median_units(rest)
    low, high = sorted((read.units, med))
    return high - low >= MIXED_RANGE_UNITS and high >= MIXED_RATIO * max(
        low, 1e-9
    )


def _lone_ok(read: ReadPoint) -> bool:
    """A single read may be an allele if it is 50+ units and unit-like."""
    return (
        read.units >= SINGLE_READ_MIN_UNITS
        and read.purity >= SINGLE_READ_MIN_PURITY
    )


def _segment(pts: Sequence[ReadPoint], notes: _Notes, depth: int = 0):
    """Splits reads into length classes by recursive accepted 2-means cuts.

    The best cut (one or more reads per side) is judged:

    * A single read far from the rest (10+ units, factor 2): a long, unit-like
      read may be an allele; any other is an outlier, set aside and named.
    * Otherwise the cut is accepted if the medians differ by at least 2 units
      and either (a) the upper median is 50+ and the lower below 50 with a
      factor of 1.5 or more (an expanded cluster next to a normal allele; no
      read-count rule), or (b) the log-space MAD test passes and, for an upper
      median below 50, the minor group holds 25% of the reads and the medians
      are ``_dprime_threshold`` pooled SDs apart.
    * A rejected cut is recorded: `CLOSE_ALLELES_UNRESOLVED` (normal-size
      groups of 2+ reads less than 5 units apart) or
      `SPLIT_NOT_SUPPORTED:<m1>/<m2>:<n_minor>`. If it failed only the 25%
      rule and the minor group is tight, that group is named
      (`LOW_FRACTION_ALLELE`) and removed instead of being merged.

    Accepted sides with a median below 50 units are segmented again; expanded
    sides are not (mosaicism makes their spread wide). Two more cuts are
    never accepted and write no flag: when both medians are 50+ units, a cut
    needs a ratio of 1.5 or more and a 25% minor group (`expanded_same`, one
    mosaic allele otherwise); inside an accepted side a cut needs 3+ reads per
    class and 5+ units (`recut`, the noise of one allele).

    Args:
      pts: Reads sorted by ascending units.
      notes: Collects outliers, named groups, and rejection flags.
      depth: 0 for the whole run, +1 inside each accepted side.

    Returns:
      A leaf (list of reads) or a pair (low node, high node).
    """
    pts = list(pts)
    if len(pts) < 2:
        return pts
    n = len(pts)
    cut = _best_cut(pts, lambda i: True)
    if cut is None:
        return pts
    low, high = pts[:cut], pts[cut:]
    lone_side = low if len(low) == 1 else high if len(high) == 1 else None
    if lone_side is not None:
        rest = high if lone_side is low else low
        if not _far_from(lone_side[0], rest):
            # An extreme read of one allele, not a candidate: try the best
            # cut with 2+ reads per side instead.
            cut = _best_cut(pts, lambda i: MIN_GROUP_READS <= i <= n - 2)
            if cut is None:
                return pts
            low, high = pts[:cut], pts[cut:]
            lone_side = None
        elif not _lone_ok(lone_side[0]):
            notes.outliers.append(lone_side[0])
            return _segment(rest, notes, depth)

    m_low, m_high = _median_units(low), _median_units(high)
    delta = m_high - m_low
    minor = low if len(low) <= len(high) else high
    reason = None
    if delta < MIN_DELTA_UNITS - _EPS:
        reason = "floor"
    elif depth >= 1 and (
        min(len(low), len(high)) < RECUT_MIN_READS or delta < CLOSE_MAX_UNITS
    ):
        # Inside a resolved allele a re-cut needs 3+ reads per class and 5+
        # units; smaller cuts there are the noise of one allele.
        reason = "recut"
    elif (
        m_high >= EXPANDED_MIN_UNITS
        and m_low < EXPANDED_MIN_UNITS
        and (m_high + 1) / (m_low + 1) >= CLEAR_RATIO
    ):
        reason = None
    elif not _accepts_log(low, high):
        reason = "log"
    elif m_low >= EXPANDED_MIN_UNITS and (
        (m_high + 1) / (m_low + 1) < CLEAR_RATIO
        or len(minor) < MIN_ALLELE_FRACTION * n
    ):
        # Both medians expanded: the spread is one mosaic allele unless the
        # clusters are far apart and balanced.
        reason = "expanded_same"
    elif m_high >= EXPANDED_MIN_UNITS:
        reason = None
    elif len(minor) < MIN_ALLELE_FRACTION * n:
        reason = "fraction"
    elif (
        delta
        < _dprime_threshold(n) * max(_pooled_sd(low, high), SD_FLOOR_UNITS)
        - _EPS
    ):
        reason = "dprime"

    if reason is None:
        sides = []
        for side in (low, high):
            if _median_units(side) < EXPANDED_MIN_UNITS and len(side) >= 2:
                sides.append(_segment(side, notes, depth + 1))
            else:
                sides.append(list(side))
        return (sides[0], sides[1])

    # Inside an accepted side only a clear candidate (5+ units) is reported:
    # smaller cuts there are the noise of a resolved allele.
    if reason in ("recut", "expanded_same"):
        return pts
    if delta > CLOSE_NOISE_UNITS and (depth == 0 or delta >= CLOSE_MAX_UNITS):
        if (
            m_high < EXPANDED_MIN_UNITS
            and min(len(low), len(high)) >= MIN_GROUP_READS
            and delta < CLOSE_MAX_UNITS
        ):
            notes.flags.append("CLOSE_ALLELES_UNRESOLVED")
        else:
            notes.flags.append(
                f"SPLIT_NOT_SUPPORTED:{m_low:.1f}/{m_high:.1f}:{len(minor)}"
            )
    if reason == "fraction" and _is_tight(minor):
        major = high if minor is low else low
        notes.named.append(list(minor))
        return _segment(major, notes, depth)
    return pts


def _leaves(node) -> list[list[ReadPoint]]:
    """All length classes under a node, lowest first."""
    if isinstance(node, tuple):
        return _leaves(node[0]) + _leaves(node[1])
    return [node]


def _segment_small(pts: Sequence[ReadPoint]):
    """Two or three reads: split at the best cut if 2+ units apart."""
    cut = _best_cut(pts, lambda i: True)
    if cut is not None:
        low, high = list(pts[:cut]), list(pts[cut:])
        if _median_units(high) - _median_units(low) >= MIN_DELTA_UNITS:
            return (low, high)
    return list(pts)


def _hp_partition(
    classes: Sequence[Sequence[ReadPoint]], pts: Sequence[ReadPoint]
):
    """Applies HP tags to reads grouped into length classes.

    Args:
      classes: Length classes of the kept reads.
      pts: Kept reads.

    Returns:
      (groups, unassigned, k, n_tagged, same_class), or None if tags cannot
      define two alleles. ``k`` counts tagged reads outside the majority
      class of their HP group; ``groups`` is meaningful only if k is 0.
      ``same_class`` is True if both HP groups have the same majority class.
    """
    hp_reads = {hp: [p for p in pts if p.hp == hp] for hp in (1, 2)}
    tagged = hp_reads[1] + hp_reads[2]
    if (
        not pts
        or len(tagged) / len(pts) < MIN_HP_FRACTION
        or min(len(hp_reads[1]), len(hp_reads[2])) < MIN_GROUP_READS
    ):
        return None
    gap = abs(_median_units(hp_reads[1]) - _median_units(hp_reads[2]))
    if (
        gap < MIN_DELTA_UNITS
        and min(len(hp_reads[1]), len(hp_reads[2])) < HP_CLOSE_MIN_READS
    ):
        return None
    class_of = {id(p): ci for ci, c in enumerate(classes) for p in c}
    majority: dict[int, int] = {}
    k = 0
    for hp, reads in hp_reads.items():
        counts = collections.Counter(class_of[id(p)] for p in reads)
        majority[hp] = sorted(counts, key=lambda c: (-counts[c], c))[0]
        k += sum(class_of[id(p)] != majority[hp] for p in reads)
    groups = {hp: list(reads) for hp, reads in hp_reads.items()}
    unassigned: list[ReadPoint] = []
    for p in pts:
        if p.hp in (1, 2):
            continue
        cls = class_of[id(p)]
        owners = [hp for hp in (1, 2) if majority[hp] == cls]
        if len(owners) == 1:
            groups[owners[0]].append(p)
        elif len(owners) == 2:
            meds = {hp: _y(_median_units(g)) for hp, g in groups.items()}
            hp = min((1, 2), key=lambda h: abs(meds[h] - _y(p.units)))
            groups[hp].append(p)
        else:
            unassigned.append(p)
    same_class = majority[1] == majority[2]
    return [groups[1], groups[2]], unassigned, k, len(tagged), same_class


def _flag_classes(classes: Sequence[Sequence[ReadPoint]]) -> str:
    return "/".join(f"{_median_units(c):.1f}x{len(c)}" for c in classes)


def _candidate_interruptions(group: Sequence[ReadPoint]) -> list[str]:
    """Reports interruption blocks seen at one offset in many reads.

    Reads with strict CTG purity below 0.5 are left out, TTG blocks are
    ignored in reads where TTG is more than half of the tokens, all blocks of
    a read are used, and positions are absolute offsets (unit positions) from
    the 3' anchor. A candidate needs the same token within 3 positions in at
    least 3 reads and at least half of the eligible reads. The result only
    annotates; it never alters a repeat count.
    """
    eligible = [p for p in group if p.strict_purity >= CANDIDATE_MIN_PURITY]
    flags = []
    if len(eligible) < CANDIDATE_MIN_READS:
        return flags
    for token in sorted(motif.INTERRUPTION_TOKENS):
        pairs = sorted(
            (block.offset, ri)
            for ri, p in enumerate(eligible)
            if not (token == "TTG" and p.ttg_fraction > CAA_MODE_TTG_FRACTION)
            for block in p.blocks
            if block.token == token
        )
        best: list[int] = []
        best_reads: set[int] = set()
        for i, (start, _) in enumerate(pairs):
            in_window = [
                (o, r)
                for o, r in pairs[i:]
                if o - start <= CANDIDATE_OFFSET_TOLERANCE
            ]
            reads = {r for _, r in in_window}
            if len(reads) > len(best_reads):
                best = [o for o, _ in in_window]
                best_reads = reads
        if (
            best
            and len(best_reads) >= CANDIDATE_MIN_READS
            and len(best_reads) >= CANDIDATE_MIN_FRACTION * len(eligible)
        ):
            flags.append(
                f"CANDIDATE_INTERRUPTION:{token}@{int(statistics.median(best))}"
            )
    return flags


def _strand_flag(group: Sequence[ReadPoint]) -> bool:
    """True if forward and reverse medians differ beyond log-space noise.

    The tolerance is the largest of 0.02 log units, a 1-unit step at the
    allele median, and three standard errors of the difference, where the
    spread is measured within each strand (MAD x 1.4826 of log values).
    """
    fwd = [p.units for p in group if not p.reverse]
    rev = [p.units for p in group if p.reverse]
    if len(fwd) < STRAND_MIN_READS or len(rev) < STRAND_MIN_READS:
        return False
    diff = abs(_y(statistics.median(fwd)) - _y(statistics.median(rev)))
    med = _median_units(group)
    one_unit = math.log((med + 2) / (med + 1))
    sd_f = MAD_SCALE * _mad([_y(u) for u in fwd])
    sd_r = MAD_SCALE * _mad([_y(u) for u in rev])
    se = math.sqrt(sd_f**2 / len(fwd) + sd_r**2 / len(rev))
    return diff > max(STRAND_MIN_LOG, one_unit, STRAND_SE_MULTIPLE * se)


def _stats(allele: int, group: Sequence[ReadPoint]) -> AlleleCall:
    """Summarizes one allele group and labels its reads.

    Args:
      allele: Allele number (1 or 2).
      group: Reads of this allele; their ``allele`` field is set.

    Returns:
      The row with median, range, MAD, percentiles (n >= 10), per-strand
      medians, remainder fraction, CAA-mode count, and per-allele flags.
    """
    values = [p.units for p in group]
    med = statistics.median(values)
    fwd = [p.units for p in group if not p.reverse]
    rev = [p.units for p in group if p.reverse]
    row = AlleleCall(
        allele=allele,
        ctg_call=math.floor(med + 0.5),
        median=med,
        minimum=min(values),
        maximum=max(values),
        n_reads=len(values),
        mad=_mad(values),
        median_fwd=statistics.median(fwd) if fwd else None,
        median_rev=statistics.median(rev) if rev else None,
        frac_rem0=sum(p.remainder == 0 for p in group) / len(group),
        motif_mode_caa=sum(
            p.ttg_fraction > CAA_MODE_TTG_FRACTION for p in group
        ),
    )
    if len(values) >= P_PERCENTILE_MIN_READS:
        row.p10, row.p90 = (
            float(x) for x in numpy.percentile(values, [10, 90])
        )
    if med >= EXPANDED_MIN_UNITS and len(values) < EXPANDED_LOW_SUPPORT_READS:
        row.flags.append("EXPANDED_LOW_SUPPORT")
    elif len(values) == 1:
        row.flags.append("SINGLE_READ_ALLELE")
    if _strand_flag(group):
        row.flags.append("STRAND_DISCORDANT")
    n_rescued = sum(p.rescued for p in group)
    if n_rescued:
        row.flags.append(f"RESCUED_READS:{n_rescued}")
    row.flags += _candidate_interruptions(group)
    for point in group:
        point.allele = allele
    return row


def call_alleles(
    points: Sequence[ReadPoint], max_lower_bound: int | None
) -> tuple[list[AlleleCall], str, list[str]]:
    """Calls up to two alleles from spanning reads.

    Reads are first cut into length classes (``_segment``). The first accepted
    cut gives the two alleles; a side that splits again into normal-size
    classes is reported as NA with `MORE_THAN_TWO_LENGTHS`. HP tags (already
    restricted to the dominant phase block by the caller) define the alleles
    when they are usable and no tagged read lies outside the majority class
    of its HP group; untagged reads follow the class of an HP group, and an
    untagged class of another length is named and left out
    (`UNTAGGED_OTHER_LENGTH`). If tagged reads disagree with the classes,
    `HP_DISCORDANT:k/n` is written and the length classes decide. If both HP
    groups fall in one class while another class exists, the tags cannot
    place that class: `HP_SAME_CLASS` is written and the lengths decide.

    Args:
      points: Spanning reads that passed QC (``allele`` is filled in).
      max_lower_bound: Longest periodic run seen in one-anchor reads (repeat
        units), or None.

    Returns:
      Two ``AlleleCall`` rows (allele 2 is NA for a single length), the
      method ("HP", "length", "none"), and run-level flags.
    """
    run_flags: list[str] = []
    n = len(points)
    if n == 0:
        rows = [AlleleCall(1), AlleleCall(2)]
        for row in rows:
            row.flags.append("NO_SPANNING_READS")
            if max_lower_bound:
                row.flags.append("LONGER_ALLELE_UNSPANNED")
        return rows, "none", run_flags
    if n < LOW_DEPTH_READS:
        run_flags.append("LOW_DEPTH")

    pts = sorted(points, key=lambda p: p.units)
    notes = _Notes()
    if n <= 3:
        tree = _segment_small(pts)
    else:
        tree = _segment(pts, notes)
    removed = {id(p) for p in notes.outliers}
    for group in notes.named:
        removed |= {id(p) for p in group}
    kept = [p for p in pts if id(p) not in removed]

    extra: list[str] = list(notes.flags)
    extra += [f"OUTLIER_READ:{p.units:g}" for p in notes.outliers]
    extra += [
        f"LOW_FRACTION_ALLELE:{_median_units(g):.1f}x{len(g)}"
        for g in notes.named
    ]
    method = "length"
    hp = _hp_partition(_leaves(tree), kept)
    sides: list = []
    if hp is not None:
        groups, unassigned, k, n_tagged, same_class = hp
        if k == 0 and same_class and len(_leaves(tree)) > 1:
            # Both haplotype groups sit in one length class while another
            # class exists (for example untagged expanded reads): the tags
            # cannot place that class, so the lengths decide.
            extra.append("HP_SAME_CLASS")
        elif k == 0:
            method = "HP"
            # HP separates what the lengths could not: not "unresolved".
            extra = [f for f in extra if f != "CLOSE_ALLELES_UNRESOLVED"]
            sides = [[g] for g in sorted(groups, key=_median_units)]
            if unassigned:
                extra.append(
                    "UNTAGGED_OTHER_LENGTH:"
                    f"{_median_units(unassigned):.1f}x{len(unassigned)}"
                )
        else:
            extra.append(f"HP_DISCORDANT:{k}/{n_tagged}")
    if not sides:
        if isinstance(tree, tuple):
            sides = [_leaves(tree[0]), _leaves(tree[1])]
        else:
            sides = [[tree]]

    rows: list[AlleleCall] = []
    for i, classes in enumerate(sides, start=1):
        if len(classes) == 1:
            rows.append(_stats(i, classes[0]))
            continue
        row = AlleleCall(i, n_reads=sum(len(c) for c in classes))
        row.flags.append(f"MORE_THAN_TWO_LENGTHS:{_flag_classes(classes)}")
        for c in classes:
            for point in c:
                point.allele = i
        rows.append(row)

    if len(rows) == 1:
        rows.append(AlleleCall(2))
        for row in rows:
            row.flags.append("SINGLE_LENGTH_OBSERVED")
        if (
            max_lower_bound is not None
            and rows[0].maximum is not None
            and max_lower_bound >= rows[0].maximum + LONGER_UNSPANNED_MARGIN
        ):
            for row in rows:
                row.flags.append("LONGER_ALLELE_UNSPANNED")
    elif (
        max_lower_bound is not None
        and rows[1].maximum is not None
        and max_lower_bound > rows[1].maximum
    ):
        rows[1].flags.append("UNSPANNED_LONGER_READS")
    for row in rows:
        row.flags += extra + run_flags
    return rows, method, run_flags
