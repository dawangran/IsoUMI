# LongUMI

LongUMI is a command-line tool for UMI correction and molecule-level deduplication on long-read single-cell BAM files.

It is designed for workflows where each read may carry:

- a cell barcode tag such as `CB`
- a raw UMI tag such as `UR`
- an optional UMI quality tag such as `UY`
- an optional gene tag such as `GX`

LongUMI groups reads by cell and transcript context, corrects low-support UMIs toward higher-support UMIs, marks duplicates, and writes a corrected BAM plus optional TSV reports.

## What It Does

Given one or more BAM files, LongUMI:

1. Splits reads into cell-barcode buckets for parallel processing.
2. Builds a molecule grouping key from:
   - cell barcode
   - gene tag, unless `--no-gene` is used
   - strand
   - splice-junction or locus structure
   - optional transcript-end bins via `--end-bin`
3. Counts UMIs inside each group.
4. Collapses low-support UMIs into high-support UMIs when they are close in sequence and satisfy count-ratio rules.
5. Optionally incorporates UMI quality information when ranking merges and computing merge confidence.
6. Writes corrected UMI tags back to BAM and optionally emits molecule, assignment, and correction-explanation tables.

## Key Features

- Long-read-aware grouping using splice-junction structure or genomic locus bins
- Multi-BAM input support
- Cell-bucket sharding for scalable parallel processing
- Quality-aware UMI scoring using `--umi-qual-tag`
- Optional merge-confidence filtering with `--min-merge-confidence`
- Optional explainability output with `--emit-explain`
- Optional transcript-end-aware grouping with `--end-bin`
- Header merging for multi-input BAMs
- Unmapped reads are passed through without participating in molecule collapse

## Build

The source code lives in [`src/`](./src). The binary is built there as `longumi`.

### Standard build

```bash
make -C src
```

### Dependencies

LongUMI depends on:

- `htslib`
- zlib / bzip2 / lzma / curl / OpenSSL libraries needed by your local `htslib`

The `Makefile` supports a few common setups:

- `pkg-config`-based `htslib`
- manual `HTSLIB_INC` and `HTSLIB_LIB`
- a fallback path for `Rhtslib` on macOS

Examples:

```bash
make -C src HTSLIB_INC=/path/to/htslib/include HTSLIB_LIB=/path/to/htslib/lib
```

If your system linker cannot find `libcrypto`, install or expose OpenSSL to the linker and rerun `make`.

## Quick Start

### Single BAM

```bash
src/longumi --bam input.bam --out sample
```

This writes:

- `sample.dedup.bam`

### Multiple BAMs

```bash
src/longumi \
  --bam sample1.bam \
  --bam sample2.bam \
  --out merged_run
```

### BAM list file

```bash
src/longumi --bam-list bam_files.txt --out run1
```

Where `bam_files.txt` contains one BAM path per line. Blank lines and `#` comments are allowed.

## Recommended Example

```bash
src/longumi \
  --bam input.bam \
  --out sample \
  --threads 8 \
  --buckets 128 \
  --emit-tsv \
  --emit-explain \
  --end-bin 50 \
  --min-merge-confidence 0.45
```

This writes:

- `sample.dedup.bam`
- `sample.molecules.tsv`
- `sample.assignments.tsv`
- `sample.corrections.tsv`

## Core Algorithm

LongUMI defines a candidate molecule using a grouping key of the form:

- `CB`
- `GX` or `GX=NA`
- reference target ID
- strand
- splice-junction hash or locus hash
- optional 5'/3' end bins when `--end-bin` is enabled

Inside each group:

1. Unique UMIs are counted.
2. UMI pairs are compared.
3. A lower-count UMI can merge into a higher-count UMI if:
   - Hamming distance is within `--ham`
   - `smaller / larger <= --ratio`
   - optional confidence threshold is satisfied
4. Connected UMIs are clustered with union-find.
5. A representative UMI is chosen per cluster.
   - higher count wins
   - if counts tie and quality-aware mode is enabled, higher mean UMI quality wins
   - if still tied, lexicographically smaller UMI wins

## Input Expectations

LongUMI expects:

- BAM inputs with compatible reference dictionaries
- UMI values stored as 2-character SAM tags
- optional UMI-quality values stored as a string tag of the same length as the UMI

For multiple BAM inputs:

- `@SQ` dictionaries must match
- metadata such as `@RG`, `@PG`, and `@CO` are merged when possible
- conflicting `@RG` or `@PG` lines with the same `ID` are treated as an error

## Important Output Semantics

### Corrected BAM

LongUMI writes a corrected UMI tag to `--umi-out` (default `UB`).

It also writes a duplicate flag tag to `--dup-flag` (default `DA`) with the following semantics:

- `DA=0`: representative read for that corrected molecule, or a passthrough read that was not deduplicated
- `DA=1`: duplicate read within the same corrected molecule

Unmapped reads:

- do not participate in molecule correction
- are passed through
- keep `UB=UR` if a raw UMI exists
- are written with `DA=0`

### Sort Order

Because reads are bucketed and concatenated, output BAM order should be treated as:

- `SO:unknown`

Do not assume coordinate-sorted or queryname-sorted output unless you sort it afterward.

## Main Parameters

### Inputs and outputs

- `--bam <FILE>`: input BAM, repeatable
- `--bam-list <FILE>`: file containing BAM paths
- `--out <PREFIX>`: output prefix

### Tags

- `--cell-tag <TAG>`: cell barcode tag, default `CB`
- `--umi-tag <TAG>`: raw UMI tag, default `UR`
- `--umi-qual-tag <TAG>`: UMI quality tag, default `UY`
- `--gene-tag <TAG>`: gene tag, default `GX`
- `--umi-out <TAG>`: corrected UMI tag to write, default `UB`
- `--dup-flag <TAG>`: duplicate flag tag, default `DA`
- `--mol-tag <TAG>`: optional molecule ID tag

### UMI correction

- `--ham <INT>`: maximum Hamming distance, default `1`
- `--ratio <FLOAT>`: require `smaller / larger <= ratio`, default `0.10`
- `--min-merge-confidence <FLOAT>`: optional confidence floor, default `0.00`
- `--no-quality-aware`: disable quality-aware ranking and quality contribution to confidence

### Grouping

- `--no-gene`: ignore gene tag during grouping
- `--locus-bin <INT>`: bin size for non-spliced reads, default `1000`
- `--sj-jitter <INT>`: splice-boundary rounding before hashing, default `10`
- `--end-bin <INT>`: add strand-aware transcript-end bins to the grouping key

### Performance

- `--threads <INT>`: worker threads, default `4`
- `--buckets <INT>`: number of cell buckets, default `64`
- `--tmp-dir <DIR>`: bucket directory, default `tmp_buckets`
- `--keep-tmp`: keep temporary bucket files

### Reports

- `--emit-tsv`: write molecule and assignment TSVs
- `--emit-explain`: write correction explanation TSV

## Report Files

### `*.molecules.tsv`

Columns:

- `key`
- `umi_corr`
- `count`
- `bucket`

This is a molecule-level summary after correction.

### `*.assignments.tsv`

Columns:

- `qname`
- `molecule_id`
- `dup`
- `bucket`

This maps each read name to a corrected molecule assignment.

### `*.corrections.tsv`

Columns:

- `key`
- `raw_umi`
- `corr_umi`
- `raw_count`
- `corr_count`
- `hamming`
- `raw_avgq`
- `corr_avgq`
- `confidence`
- `reason`
- `bucket`

`reason` is currently one of:

- `self`: no correction was needed
- `count+distance`: correction supported by count and sequence distance
- `count+quality`: correction additionally supported by quality information

## Practical Tips

- Start with defaults if your BAM already carries `CB`, `UR`, `UY`, and `GX`.
- Use `--no-gene` only if gene tags are absent or unreliable.
- Use `--end-bin` when you want molecule grouping to better respect transcript ends.
- Use `--emit-explain` when benchmarking or tuning parameters.
- Raise `--buckets` for larger datasets to spread memory across more buckets.
- Tighten `--min-merge-confidence` if you want more conservative UMI correction.

## Example Workflows

### Conservative correction

```bash
src/longumi \
  --bam input.bam \
  --out conservative \
  --ham 1 \
  --ratio 0.05 \
  --min-merge-confidence 0.60 \
  --emit-explain
```

### Isoform-aware grouping

```bash
src/longumi \
  --bam input.bam \
  --out isoform_run \
  --end-bin 25 \
  --emit-tsv \
  --emit-explain
```

### Tag remapping

```bash
src/longumi \
  --bam input.bam \
  --out custom_tags \
  --cell-tag XC \
  --umi-tag XM \
  --umi-qual-tag XQ \
  --gene-tag GN \
  --umi-out UB \
  --dup-flag DA
```

## Current Limitations

- UMI distance is still sequence-based and currently uses Hamming distance only.
- UMI quality is used for ranking and confidence scoring, not a full probabilistic error model.
- Final BAM linking may depend on your local OpenSSL and `htslib` setup.
- Output BAM is not coordinate-sorted by construction.

## Repository Layout

```text
.
├── README.md
└── src/
    ├── Makefile
    ├── main.c
    ├── cli.c
    ├── pipeline.c
    ├── key_build.c
    ├── sj.c
    ├── umi.c
    └── ...
```

## Version

Current version string in the source:

- `LongUMI 0.3.2 (fixed pipeline.c)`
