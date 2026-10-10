# dmpkCTGcounter

dmpkCTGcounter counts the CTG repeat units of both alleles at the (CTG)n locus in the 3' untranslated region of *DMPK* (GRCh38) from one BAM or CRAM alignment file of Oxford Nanopore Technologies (ONT) reads. The input can come from adaptive sampling or from standard whole-genome sequencing. The tool does one thing: it reports a median repeat count, a range, and a read count for each allele, with flags that state how far the call can be trusted. It does not call variants, methylation, or interruptions, and it gives no clinical interpretation. The tool is for research use only. This tool is an in-house script provided as supplementary material, with no claim of validation or support.

* Version: 0.1.0 (Python 3.12)
* Reference build: GRCh38 only; other builds stop the run with exit code 2.
* Locus: NC_000019.10:g.45770205_45770264, NM_004409.5:c.\*224_\*283CTG[20]. *DMPK* lies on the minus strand, so the GRCh38 plus strand carries the tract as (CAG)n and the gene strand as (CTG)n. Counts are reported as CTG units.

## Method

Each read is measured from its own sequence, never from its alignment CIGAR.

1. Read selection: One fetch of chr19:45,760,000-45,780,000 (1-based). Records with flags 0x4, 0x100, 0x200, or 0x400 and duplex parent reads (`dx` tag -1) are dropped. The primary SEQ is searched whole, soft-clipped bases included. A supplementary record is resolved to its primary through the SA tag, and the primary is reverse-complemented when the two records lie on opposite strands. A primary with hard clips is skipped and counted.
2. Anchors: Two fixed 40 bp anchors on the plus strand flank the tract: chr19:45,770,160-45,770,199 and chr19:45,770,267-45,770,306 (GRCh38, 1-based). Both are located with edlib (at most 6 edits). A second hit within max(6, best + 3) edits rejects the read as ambiguous. The reverse complement is searched only to reject inverted chimeras and reverse-complemented records.
3. Repeat units: In the reference the gap between the anchors is 67 bp: a 5 bp CCCCC margin, the 60 bp tract (CAG)20, and a 2 bp CA margin. Units = (gap - 7) / 3, a real number, which gives 20.0 for the reference. Interruptions (CCG, CGG, CTC, CAG, TTG) are counted as length. The remainder (gap - 7) mod 3 is recorded and never used as a filter.
4. Anchor rescue: If one anchor is found and its periodic run ends at least 30 bp before the read end, the other anchor is searched as a 100 bp anchor. The 100 bp anchors are chr19:45,770,100-45,770,199 (left) and chr19:45,770,267-45,770,366 (right). The search allows 25 edits and covers 50 bp before to 250 bp after the end of the run. Rescued reads are counted and flagged `RESCUED_READS:k`.
5. One-anchor reads: A read with one anchor whose periodic run reaches the read end gives a lower bound on the repeat count. These bounds are reported (`n_one_anchor`, `max_lower_bound`) and never enter an allele call.
6. Purity filter: A read with at least 20 units is rejected when fewer than 50% of its tract triplets are CTG or a one-base variant (`--min-purity`, `--purity-min-units`).
7. **Allele calling.**
    * Reads are cut into length classes in log(units + 1) space by repeated 2-means cuts. A cut is accepted when the medians differ by at least 2 units and one of two conditions holds. Either the upper median is at least 50 units and at least 1.5 times a lower median below 50, or the log-space difference is at least three pooled SDs estimated from the median absolute deviation (MAD x 1.4826). For two normal-size groups (upper median below 50) the smaller group must also hold at least 25% of the reads. The medians must differ by at least 6.0, 5.5, 5.0, or 4.5 pooled SDs for fewer than 12, fewer than 20, fewer than 40, or 40 or more reads. These rules keep a noisy single allele from being split.
    * The first accepted cut gives the two alleles. If a side below 50 units splits again, that allele is NA and carries `MORE_THAN_TWO_LENGTHS`; a median is never reported across accepted classes.
    * One length class is never called homozygous: allele 2 is NA with `SINGLE_LENGTH_OBSERVED`.
    * HP (haplotype) tags are used only within the dominant PS (phase set) block among the spanning reads. They define the alleles when three conditions hold: at least 60% of reads are tagged, HP 1 and HP 2 each have at least 2 reads (3 when their medians are less than 2 units apart), and every tagged read lies in the majority length class of its HP group. Otherwise the length classes decide.
    * Each allele row reports the median, minimum, and maximum of its spanning reads, `ctg_call` = floor(median + 0.5), and the read count.

The procedure is deterministic and uses no random numbers. Thresholds are in reads, units, and log ratios and do not assume a coverage level.

## Installation

The tool runs from an Apptainer image built from `env/dmpk.def`. No conda or pip installation is needed on the host. The image is not stored in the repository; build it once from the repository root:

```bash
apptainer build env/dmpk.sif env/dmpk.def
```

The definition pins every component:

| Component | Version |
|---|---|
| Base image | python:3.12.12-slim-bookworm (pinned by digest) |
| pip | 26.1.2 |
| pysam | 0.23.3 |
| edlib | 1.3.9.post1 |
| numpy | 2.4.6 |

`env/requirements.txt` lists the exact pins. The build runs `pip check` and writes `pip freeze` to `/opt/provenance/pip-freeze.txt` inside the image.

## Usage

```bash
apptainer exec --bind /path/to/data env/dmpk.sif \
    dmpkCTGcounter sample.cram --reference GRCh38.fasta --out sample.ctg.tsv
```

`apptainer run env/dmpk.sif <arguments>` is equivalent. Inside the image the entry point is `dmpkCTGcounter`, or `python -m dmpkctgcounter`.

```text
usage: dmpkCTGcounter [-h] [--reference REFERENCE] [--out OUT]
                      [--sample SAMPLE] [--basecall-model BASECALL_MODEL]
                      [--min-mapq MIN_MAPQ] [--max-edits MAX_EDITS]
                      [--min-purity MIN_PURITY]
                      [--purity-min-units PURITY_MIN_UNITS] [--per-read TSV]
                      [--size-class] [--version]
                      alignment
```

| Option | Meaning |
|---|---|
| `alignment` | One BAM or CRAM with an index (.bai, .csi, or .crai), GRCh38, ONT. |
| `--reference` | Indexed GRCh38 FASTA (.fai required). It verifies the anchor constants, the margins, and the (CAG)20 tract, and it decodes CRAM. Default: the FASTA named by `DMPKCTGCOUNTER_REFERENCE`, if set. |
| `--out` | Output TSV (default: stdout). A one-line summary and any warnings go to stderr. |
| `--sample` | Sample name (default: `@RG SM`, or the file name without extension). |
| `--basecall-model` | Basecall model name (default: read from the `@PG` and `@RG DS` header lines). |
| `--min-mapq` | Skip primary records below this MAPQ (default 0). |
| `--max-edits` | Edit budget per 40 bp anchor (default 6). |
| `--min-purity` | Minimum unit-like token fraction (default 0.5). |
| `--purity-min-units` | The purity filter applies from this many units (default 20). |
| `--per-read` | Also write one row per read to this TSV. |
| `--size-class` | Add the `size_class_research_use` column. |
| `--version` | Print the version. |

Reference handling:

* A CRAM file needs the reference it was encoded against. Without `--reference` the tool stops with exit code 2 if the CRAM cannot be decoded.
* For a BAM file the reference is optional. Without it the built-in anchor constants are used, and the TSV header states that no reference was checked.
* If `--reference` is omitted and the environment variable `DMPKCTGCOUNTER_REFERENCE` names an indexed FASTA, that file is used and the TSV header says "from DMPKCTGCOUNTER_REFERENCE". If the variable is set but the `.fai` index is missing, the tool stops with exit code 2. If neither is given, the built-in anchor constants apply. Pass `--reference` explicitly for reproducible runs.
* The `@SQ M5` of chr19 is compared with the FASTA; a mismatch gives a warning.
* The tool disables the `REF_PATH` and network reference lookups of htslib. A local `@SQ UR` path in a CRAM header is still used by htslib if that file exists.

Exit code 0 on success, including runs whose rows are all NA. Exit code 2 for input errors: wrong genome build (chr19 length not 58,617,616, or chr1 length not 248,956,422 when present), missing index, missing or mismatched reference, undecodable records, or an unreadable file.

## Output

The TSV starts with a `#` comment block: tool version, command line, package versions, input and reference status, build, tract, anchors, offset, basecall model, parameters, method, read counts by status and rejection reason, and warnings. A column line and two rows (allele 1, the shorter, and allele 2, the longer) follow.

Illustrative synthetic output (values are invented, comment block and some columns omitted):

```text
sample   allele  ctg_call  ctg_median  ctg_min  ctg_max  n_reads  flags
example  1       13        13.0        12.3     13.7     14       .
example  2       612       612.3       540.0    701.7    4        EXPANDED_LOW_SUPPORT
```

| Column | Meaning |
|---|---|
| `sample` | `@RG SM`, the file stem, or `--sample`. |
| `allele` | 1 (shorter) or 2 (longer). |
| `ctg_call` | floor(median + 0.5); NA if not called. |
| `ctg_median`, `ctg_min`, `ctg_max` | Median, minimum, and maximum of the real-valued units of the spanning reads. |
| `n_reads` | Spanning reads assigned to the allele. |
| `flags` | Semicolon-separated flags; `.` if none. |
| `mad` | Median absolute deviation of the units. |
| `p10`, `p90` | 10th and 90th percentiles; NA below 10 reads. |
| `median_fwd`, `median_rev` | Median units of forward-strand and reverse-strand reads. |
| `n_one_anchor`, `max_lower_bound` | One-anchor reads of the run and the longest periodic run among them (units). Same value on both rows. |
| `frac_rem0` | Fraction of the allele's reads with (gap - 7) divisible by 3. |
| `motif_mode_caa` | Reads of the allele in CAA mode (TTG tokens above 50% of unit-like tokens). |
| `size_class_research_use` | With `--size-class` only; see below. |

### Flags

Flags that describe the whole run are written on both rows.

| Flag | Meaning |
|---|---|
| `NO_READS_IN_REGION` | No record in chr19:45,760,000-45,780,000. |
| `NO_USABLE_READS` | Records present, but all were skipped or excluded before measurement (flags, duplex parent, hard clips, MAPQ). |
| `NO_SPANNING_READS` | Reads were measured; none span the tract. |
| `LOW_DEPTH` | Fewer than 5 spanning reads. |
| `SINGLE_LENGTH_OBSERVED` | Only one length class; allele 2 is NA. |
| `SINGLE_READ_ALLELE` | The allele rests on one read below 50 units. |
| `EXPANDED_LOW_SUPPORT` | Median of 50 units or more from fewer than 5 reads. |
| `OUTLIER_READ:<units>` | A single read at least 10 units and a factor of 2 from the rest, below 50 units or with a unit-like fraction below 0.8; left out of the statistics. |
| `SPLIT_NOT_SUPPORTED:<m1>/<m2>:<n>` | The best cut was rejected; medians of both groups and reads in the smaller group. |
| `CLOSE_ALLELES_UNRESOLVED` | Normal-size groups of 2 or more reads less than 5 units apart were not separated. Not written when HP tags resolve them. |
| `LOW_FRACTION_ALLELE:<median>x<n>` | A tight group below the 25% rule; not an allele, no size class, not merged into the main allele. |
| `MORE_THAN_TWO_LENGTHS:<median>x<n>/...` | A side split into two or more normal-size classes; the allele is NA and each class is listed. |
| `UNTAGGED_OTHER_LENGTH:<median>x<n>` | Untagged reads of a length that matches neither HP group; left out. |
| `HP_DISCORDANT:k/n` | k of n tagged reads lie outside the majority length class of their HP group; length classes decide. |
| `HP_SAME_CLASS` | Both HP groups sit in one length class while another class exists (for example untagged expanded reads). The tags cannot place that class, so they are ignored and lengths decide. |
| `STRAND_DISCORDANT` | Forward and reverse medians differ beyond a log-space tolerance (at least 3 reads per strand). |
| `RESCUED_READS:k` | k reads of the allele used a rescued anchor. |
| `CANDIDATE_INTERRUPTION:<token>@<offset>` | The same interruption block (CCG, CGG, CTC, CAG, or TTG) at the same offset (within 3 positions) from the 3' anchor in at least 3 reads and half of the eligible reads. Annotation only; it never changes a count. |
| `LONGER_ALLELE_UNSPANNED` | One-anchor reads show a run at least 5 units longer than the single observed allele, or exist when no read spans. |
| `UNSPANNED_LONGER_READS` | A one-anchor run is longer than the maximum of allele 2. |

### Size class

With `--size-class` the tool adds `size_class_research_use`: `<5`, `normal(5-34)`, `variable(35-49)`, or `fully_variable(>=50)`, assigned from `ctg_call`. The thresholds follow the 2024 ACMG technical standard for myotonic dystrophy type 1 testing (Seifert et al., 2024). The value is NA when there is no call or the allele rests on one read below 50 units. The column is a research label and not a diagnostic classification.

## License

MIT; see [LICENSE](LICENSE).

## Citation

The software supports the *DMPK* repeat counts in an unpublished manuscript. Please cite:

Tsai, M.-J. M., et al. (n.d.). *Congenital myotonic dystrophy type 1 as a sentinel diagnosis: Intergenerational epigenetic divergence and underrecognized maternal DM1* [Manuscript in preparation]. National Taiwan University Hospital.

A full reference will replace this entry upon publication. 

Software maintained by Jen-Hsiang Ou [oujenhsiang@ntu.edu.tw], Department of Medical Genetics, National Taiwan University Hospital, maintains the software (version 0.1.0).

## References

* Seifert, B. A., Reddi, H. V., Kang, B. E., Bean, L. J. H., Shealy, A., & Rose, N. C. (2024). Myotonic dystrophy type 1 testing, 2024 revision: A technical standard of the American College of Medical Genetics and Genomics (ACMG). *Genetics in Medicine, 26*(8), 101145. https://doi.org/10.1016/j.gim.2024.101145
