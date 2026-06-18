# IsoUMI Parameter Guide

This guide explains how each IsoUMI parameter affects grouping, UMI correction,
deduplication, and reproducibility. The short help remains available with
`src/isoumi --help`.

## Inputs And Outputs

### `--bam <FILE>`

Input SAM/BAM file. This option is repeatable.

Use repeated `--bam` arguments when one run should combine multiple input files
before molecule correction:

```bash
src/isoumi --bam sample1.bam --bam sample2.bam --out merged
```

For multiple inputs, reference dictionaries must match exactly. Compatible
metadata lines are merged; conflicting `@RG` or `@PG` IDs are rejected.

By default, repeated inputs are deduplicated together. This is appropriate for
technical splits or lanes from the same cell-barcode namespace. For independent
samples or libraries that may reuse cell barcodes, use `--isolate-inputs` or
include an existing source boundary with `--source-tag`.

### `--bam-list <FILE>`

Text file containing one SAM/BAM path per line. Blank lines and lines beginning
with `#` are ignored.

This is useful for large batches where repeated `--bam` arguments become hard
to read.

### `--out <PREFIX>`

Required output prefix. The corrected BAM is written as:

```text
<PREFIX>.dedup.bam
```

Optional reports use the same prefix:

```text
<PREFIX>.molecules.tsv
<PREFIX>.assignments.tsv
<PREFIX>.corrections.tsv
```

## Tag Parameters

All tag-name parameters must be two-character SAM tags.

### `--cell-tag <TAG>` default `CB`

Cell barcode tag. Reads without this tag are passed through and do not
participate in UMI correction, molecule reports, or assignment reports.

Changing this is necessary for non-10x or custom preprocessing pipelines.

### `--umi-tag <TAG>` default `UR`

Raw UMI tag read from each input record. Reads without this tag are passed
through with duplicate flag `0`.

UMI correction compares strings with Hamming distance, so UMIs must have equal
length to be compared.

### `--umi-qual-tag <TAG>` default `UY`

Optional UMI quality string. The string should have the same length as the UMI
and use Phred+33 characters.

When available, UMI quality contributes to seed ranking and merge-confidence
scoring. Missing or length-mismatched quality strings are ignored for that UMI.

### `--gene-tag <TAG>` default `GX`

Gene tag used in the grouping key unless `--no-gene` is set. Reads missing the
gene tag are grouped with `GX=NA`.

Use a gene tag that matches the annotation strategy used upstream.

### `--source-tag <TAG>`

Optional existing SAM tag to include in the grouping key. A common choice is
`RG` when read groups encode sample, library, or lane boundaries that should not
be deduplicated together.

Reads missing this tag are grouped with `SRC=NA`.

### `--input-scope-tag <TAG>` default `zi`

Temporary tag used internally by `--isolate-inputs`. The tag is added while
writing bucket BAMs, included in the grouping key as `IN=...`, and removed
before the final BAM is written.

The tag must not already exist in input records. If it does, choose a different
two-character tag.

### `--umi-out <TAG>` default `UB`

Output tag for the corrected UMI. For reads that are not corrected but have a
raw UMI, this tag is set to the raw UMI.

### `--dup-flag <TAG>` default `DA`

Integer duplicate flag written to each output record:

- `0`: representative read for a corrected molecule, or passthrough read
- `1`: duplicate read within a corrected molecule

Representative reads are selected by mapping quality, aligned reference span,
query length, and then stable input order.

IsoUMI also sets or clears the standard SAM duplicate bit (`0x400`) by default
so downstream BAM tools can recognize duplicates without reading `DA`.

### `--mol-tag <TAG>`

Optional output tag containing a molecule identifier:

```text
CB|grouping_key|corrected_UMI
```

This is useful when inspecting BAM records directly. The same identifier appears
in `*.assignments.tsv` when `--emit-tsv` is enabled.

## Grouping Parameters

Grouping parameters decide which reads are allowed to compare UMIs. They do not
directly change UMI distance calculations.

### `--no-gene`

Ignore the gene tag during grouping. The grouping key uses `GX=NA` for all
reads.

Use this only when gene tags are absent, unreliable, or intentionally excluded
from an analysis. It makes groups broader and can increase UMI comparisons.

### `--no-structure`

Baseline mode for application-note comparisons. It ignores splice-junction,
locus, and transcript-end structure in the grouping key.

With this option, grouping still uses:

```text
cell barcode + gene/NA + optional source/input scope + reference target + strand
```

This mode is intentionally less isoform-aware and is useful for showing how much
splice-junction and locus context prevent over-collapse.

### `--isolate-inputs`

Treat each `--bam` argument or `--bam-list` entry as a separate source during
grouping. IsoUMI writes an internal input-scope tag into temporary bucket BAMs,
adds it to the molecule key, and removes it before writing the final BAM.

Use this when separate samples or libraries may reuse the same cell barcode
values. Do not use it for technical splits that should be deduplicated together
across files.

### `--locus-bin <INT>` default `1000`

Bin size in base pairs for non-spliced reads, defined as reads without a CIGAR
`N` operation.

For non-spliced reads, IsoUMI bins the reference start and end positions:

```text
start_bin = floor(start / locus_bin) * locus_bin
end_bin   = floor(end   / locus_bin) * locus_bin
```

The pair is hashed into an `LX=...` grouping component.

Smaller values are more precise and reduce over-grouping. Larger values tolerate
more alignment-end variation but can merge nearby contexts.

Typical choices:

- `1000`: default broad locus grouping for noisy long-read alignments
- `100` or `250`: stricter locus separation
- `2000`: more permissive grouping for sparse data

### `--sj-jitter <INT>` default `10`

Rounding window for splice-junction boundaries before hashing. It applies only
to spliced reads with CIGAR `N` operations.

Each junction start/end is rounded to the nearest multiple of `sj-jitter`.
For example, with `--sj-jitter 10`:

```text
1004 -> 1000
1005 -> 1010
1009 -> 1010
```

Lower values are stricter. Higher values tolerate small alignment boundary
differences but can merge nearby alternative splice sites.

Typical choices:

- `0` or `1`: exact junction boundaries
- `10`: default tolerance for small long-read alignment jitter
- `20` or `50`: more permissive; inspect results carefully

### `--end-bin <INT>`

Adds strand-aware transcript 5' and 3' end bins to the grouping key.

For positive-strand reads:

```text
E5 = read_start
E3 = read_end
```

For negative-strand reads:

```text
E5 = read_end
E3 = read_start
```

Both ends are binned with `floor(position / end_bin) * end_bin`.

Use this when transcript-end differences are biologically important. Smaller
values are more isoform-specific but more sensitive to truncation, soft
clipping, and alignment variation.

Typical choices:

- disabled: default, groups by junction or locus structure only
- `25` or `50`: transcript-end-aware grouping
- `100` or `200`: softer end-aware grouping for noisy data

## UMI Correction Parameters

UMI correction is performed inside each grouping key. Candidate seed UMIs are
processed from strongest to weakest. A raw UMI is corrected to a seed only when
that raw-to-seed pair directly satisfies all enabled filters.

### `--ham <INT>` default `1`

Maximum Hamming distance between the raw UMI and corrected seed UMI.

`--ham 1` means the final reported `raw_umi -> corr_umi` correction must differ
by at most one position. Chain-based corrections that exceed this direct
distance are not allowed.

Use `0` to disable sequence-error correction while still marking exact-UMI
duplicates.

### `--ratio <FLOAT>` default `0.10`, range `0..1`

Maximum allowed support ratio:

```text
raw_count / seed_count <= ratio
```

Lower values are more conservative. Higher values correct more aggressively and
can merge true low-abundance molecules if the grouping context is broad.

Typical choices:

- `0.05`: conservative
- `0.10`: default
- `0.50` or `0.60`: useful for small tests or highly error-prone UMIs
- `1.00`: permits equal-count correction according to seed ranking

### `--min-merge-confidence <FLOAT>` default `0.00`, range `0..1`

Optional heuristic confidence floor for corrections. The score combines count
ratio, UMI distance, and optional quality evidence. It is not a posterior
probability.

Set this above zero when you want a stricter correction policy:

```bash
--min-merge-confidence 0.45
--min-merge-confidence 0.60
```

### `--no-quality-aware`

Disable quality-aware seed ranking and quality contribution to merge
confidence.

Use this when UMI quality tags are absent, unreliable, or not comparable across
inputs.

## Performance And Temporary Files

### `--threads <INT>` default `4`

Number of worker threads for per-bucket processing. Threading happens across
cell-barcode buckets.

If IsoUMI is built without OpenMP support, this option is accepted but
per-bucket processing runs serially. The program prints a runtime notice in
that case.

### `--buckets <INT>` default `64`

Number of cell-barcode buckets. Increasing this can lower peak memory per bucket
on large datasets, but creates more temporary files.

Typical choices:

- `16` or `32`: small datasets
- `64`: default
- `128` or `256`: larger datasets or high-depth cells

### `--tmp-dir <DIR>` default `<out>.isoumi.tmp.<pid>`

Directory for bucket BAMs and temporary per-bucket reports. The directory must
not already exist; this prevents accidental reuse of stale bucket files.

### `--keep-tmp`

Keep temporary bucket files after the run. Use this for debugging bucket-level
failures or inspecting intermediate BAMs.

## Report Parameters

### `--emit-tsv`

Write molecule and assignment reports:

```text
<out>.molecules.tsv
<out>.assignments.tsv
```

Use this for benchmarking, molecule count comparisons, and debugging duplicate
labels.

### `--emit-explain`

Write correction explanations:

```text
<out>.corrections.tsv
```

This is the most useful report for parameter tuning. It shows raw UMI,
corrected UMI, counts, Hamming distance, quality summaries, confidence, reason,
and bucket.

The correction report uses `seed_count` and `seed_avgq` for the target UMI
before lower-support UMIs are merged into it. Final molecule read counts are in
`<out>.molecules.tsv`.

### `--no-bam-dup-flag`

Leave the standard SAM duplicate bit unchanged. The custom duplicate tag
selected by `--dup-flag` is still written.

## Recommended Profiles

### Default exploratory run

```bash
src/isoumi --bam input.bam --out sample --emit-tsv --emit-explain
```

### More conservative correction

```bash
src/isoumi \
  --bam input.bam \
  --out conservative \
  --ham 1 \
  --ratio 0.05 \
  --min-merge-confidence 0.60 \
  --emit-explain
```

### Transcript-end-aware grouping

```bash
src/isoumi \
  --bam input.bam \
  --out end_aware \
  --end-bin 50 \
  --emit-tsv \
  --emit-explain
```

### Application-note baseline

```bash
src/isoumi \
  --bam input.bam \
  --out baseline_no_structure \
  --no-structure \
  --emit-tsv \
  --emit-explain
```
