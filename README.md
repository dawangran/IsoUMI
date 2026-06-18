# IsoUMI

IsoUMI is a command-line tool for UMI correction and molecule-level deduplication on long-read single-cell BAM files.

It is designed for workflows where each read may carry:

- a cell barcode tag such as `CB`
- a raw UMI tag such as `UR`
- an optional UMI quality tag such as `UY`
- an optional gene tag such as `GX`

IsoUMI groups reads by cell and transcript context, corrects low-support UMIs toward higher-support UMIs, marks duplicates, and writes a corrected BAM plus optional TSV reports.

## What It Does

Given one or more BAM files, IsoUMI:

1. Splits reads into cell-barcode buckets for parallel processing.
2. Builds a molecule grouping key from:
   - cell barcode
   - gene tag, unless `--no-gene` is used
   - optional source tag or input scope
   - strand
   - splice-junction or locus structure, unless `--no-structure` is used
   - optional transcript-end bins via `--end-bin`
3. Counts UMIs inside each group.
4. Collapses low-support UMIs directly into higher-support seed UMIs when they are close in sequence and satisfy count-ratio rules.
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
- Baseline no-structure grouping for application-note comparisons via `--no-structure`
- Optional source-aware grouping via `--source-tag` or `--isolate-inputs`
- Standard SAM duplicate flag support, with custom duplicate tag retained
- Header merging for multi-input BAMs
- Unmapped reads are passed through without participating in molecule collapse

## Availability And Citation

IsoUMI is released under the MIT license. The publication release is version `0.1.0`.

- Source code: update `CITATION.cff` with the final GitHub URL before submission.
- Archive DOI: add the Zenodo or institutional archive DOI after tagging the release.
- Citation metadata: see [`CITATION.cff`](./CITATION.cff).
- Application note draft: see [`docs/application_note.md`](./docs/application_note.md).
- Manuscript and figure plan: see [`docs/application_note_manuscript_plan.md`](./docs/application_note_manuscript_plan.md) and [`docs/figure_design.md`](./docs/figure_design.md).

## Build

The source code lives in [`src/`](./src). The binary is built there as `isoumi`.

### Standard build

```bash
make
```

The default build avoids CPU-specific `-march=native` flags so release builds
remain portable across machines. For a local performance build on the same CPU
family, use:

```bash
make NATIVE=1
```

### Dependencies

IsoUMI depends on:

- `htslib`
- zlib / bzip2 / lzma / curl / OpenSSL libraries needed by your local `htslib`

The `Makefile` supports a few common setups:

- `pkg-config`-based `htslib`
- manual `HTSLIB_INC` and `HTSLIB_LIB`
- a fallback path for `Rhtslib` on macOS

On macOS, OpenMP is not enabled by default because the system compiler commonly lacks it. In that case `--threads` is accepted but bucket processing runs serially and IsoUMI prints a runtime notice. Build with appropriate `OPENMP_CFLAGS` and `OPENMP_LDFLAGS` if you install an OpenMP-capable toolchain.

Examples:

```bash
make -C src HTSLIB_INC=/path/to/htslib/include HTSLIB_LIB=/path/to/htslib/lib
```

If your `htslib` build needs OpenSSL symbols, pass them explicitly through `CRYPTO_LIBS`, for example:

```bash
make CRYPTO_LIBS="-lcrypto"
```

Run the smoke test after building:

```bash
make test
```

Continuous integration is configured in [`.github/workflows/ci.yml`](./.github/workflows/ci.yml) for Ubuntu builds with `libhts-dev`.

## Quick Start

### Single BAM

```bash
src/isoumi --bam input.bam --out sample
```

This writes:

- `sample.dedup.bam`

### Multiple BAMs

```bash
src/isoumi \
  --bam sample1.bam \
  --bam sample2.bam \
  --out merged_run
```

### BAM list file

```bash
src/isoumi --bam-list bam_files.txt --out run1
```

Where `bam_files.txt` contains one BAM path per line. Blank lines and `#` comments are allowed.

### Minimal Reproducible Example

A tiny SAM fixture is provided under [`examples/`](./examples) for checking command-line behavior:

```bash
make
sh examples/run_minimal.sh
```

The example writes outputs under `examples/output/` and should show raw UMI `AAAT` collapsing to `AAAA` in `minimal.corrections.tsv`.

## Recommended Example

```bash
src/isoumi \
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

IsoUMI defines a candidate molecule using a grouping key of the form:

- `CB`
- `GX` or `GX=NA`
- optional `SRC` from `--source-tag`
- optional `IN` from `--isolate-inputs`
- reference target ID
- strand
- splice-junction hash or locus hash
- optional 5'/3' end bins when `--end-bin` is enabled

When `--no-structure` is used, the key intentionally omits splice-junction, locus, and transcript-end terms. This mode is intended as an internal baseline for evaluating the value of isoform-aware grouping.

Inside each group:

1. Unique UMIs are counted.
2. Candidate seed/raw UMI pairs are compared.
3. A lower-count raw UMI can merge into a higher-count seed UMI if:
   - Hamming distance is within `--ham`
   - `smaller / larger <= --ratio`
   - optional confidence threshold is satisfied
4. Candidate seed UMIs are processed from strongest to weakest.
5. A raw UMI is corrected only to a seed UMI that directly satisfies the distance, count-ratio, and confidence filters.
6. Seed priority is:
   - higher count wins
   - if counts tie and quality-aware mode is enabled, higher mean UMI quality wins
   - if still tied, lexicographically smaller UMI wins

## Input Expectations

IsoUMI expects:

- BAM inputs with compatible reference dictionaries
- UMI values stored as 2-character SAM tags
- cell barcode values stored in `--cell-tag`
- optional UMI-quality values stored as a string tag of the same length as the UMI

For multiple BAM inputs:

- `@SQ` dictionaries must match
- metadata such as `@RG`, `@PG`, and `@CO` are merged when possible
- conflicting `@RG` or `@PG` lines with the same `ID` are treated as an error
- inputs are deduplicated together by default; use `--isolate-inputs` for independent samples or `--source-tag RG` when read-group tags define the desired source/library boundary

## Important Output Semantics

### Corrected BAM

IsoUMI writes a corrected UMI tag to `--umi-out` (default `UB`).

It also writes a duplicate flag tag to `--dup-flag` (default `DA`) with the following semantics:

- `DA=0`: best representative read for that corrected molecule, or a passthrough read that was not deduplicated
- `DA=1`: duplicate read within the same corrected molecule

By default, IsoUMI also sets the standard SAM duplicate bit (`0x400`) consistently with `DA`. Use `--no-bam-dup-flag` if you want to leave the standard BAM flag untouched and rely only on the custom tag.

For deduplicated molecules, the representative read is chosen by higher mapping quality, then longer aligned reference span, then longer query length, then earlier input order as a stable tie-breaker.

Unmapped reads:

- do not participate in molecule correction
- are passed through
- keep `UB=UR` if a raw UMI exists
- are written with `DA=0`

Mapped reads without a cell barcode:

- do not participate in molecule correction
- are passed through
- keep `UB=UR` if a raw UMI exists
- are written with `DA=0`
- are omitted from molecule, assignment, and correction reports

### Sort Order

Because reads are bucketed and concatenated, output BAM order should be treated as:

- `SO:unknown`

Do not assume coordinate-sorted or queryname-sorted output unless you sort it afterward.

## Main Parameters

This section is a quick reference. For detailed parameter logic, tuning
recommendations, and interactions, see [`docs/parameters.md`](./docs/parameters.md).

### Inputs and outputs

- `--bam <FILE>`: input BAM, repeatable
- `--bam-list <FILE>`: file containing BAM paths
- `--out <PREFIX>`: output prefix

### Tags

- `--cell-tag <TAG>`: cell barcode tag, default `CB`
- `--umi-tag <TAG>`: raw UMI tag, default `UR`
- `--umi-qual-tag <TAG>`: UMI quality tag, default `UY`
- `--gene-tag <TAG>`: gene tag, default `GX`
- `--source-tag <TAG>`: existing source/library tag to include in grouping, for example `RG`
- `--umi-out <TAG>`: corrected UMI tag to write, default `UB`
- `--dup-flag <TAG>`: duplicate flag tag, default `DA`
- `--mol-tag <TAG>`: optional molecule ID tag
- `--input-scope-tag <TAG>`: temporary internal tag used by `--isolate-inputs`, default `zi`

### UMI correction

- `--ham <INT>`: maximum Hamming distance, default `1`
- `--ratio <FLOAT>`: require `smaller / larger <= ratio`, range `0..1`, default `0.10`
- `--min-merge-confidence <FLOAT>`: optional confidence floor, default `0.00`
- `--no-quality-aware`: disable quality-aware ranking and quality contribution to confidence

### Grouping

- `--no-gene`: ignore gene tag during grouping
- `--no-structure`: ignore splice-junction, locus, and transcript-end structure during grouping; useful as a baseline comparison mode
- `--isolate-inputs`: treat each `--bam` or `--bam-list` entry as a separate grouping source
- `--locus-bin <INT>`: bin size for non-spliced reads, default `1000`
- `--sj-jitter <INT>`: splice-boundary rounding to the nearest jitter multiple before hashing, default `10`
- `--end-bin <INT>`: add strand-aware transcript-end bins to the grouping key

### Performance

- `--threads <INT>`: worker threads, default `4`
- `--buckets <INT>`: number of cell buckets, default `64`
- `--tmp-dir <DIR>`: bucket directory, default `<out>.isoumi.tmp.<pid>`. The directory must not already exist.
- `--keep-tmp`: keep temporary bucket files

### Reports

- `--emit-tsv`: write molecule and assignment TSVs
- `--emit-explain`: write correction explanation TSV
- `--no-bam-dup-flag`: do not update the standard SAM duplicate bit

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
- `seed_count`
- `hamming`
- `raw_avgq`
- `seed_avgq`
- `confidence`
- `reason`
- `bucket`

`reason` is currently one of:

- `self`: no correction was needed
- `count+distance`: correction supported by count and sequence distance
- `count+quality`: correction additionally supported by quality information

`seed_count` and `seed_avgq` describe the target seed UMI before lower-support UMIs are merged into it. The final corrected molecule read count is reported in `*.molecules.tsv`.

For `self` rows, `confidence` is reported as `1.0`. For actual UMI corrections, `confidence` is a heuristic merge score derived from count ratio, Hamming distance, and optional UMI quality; it is not a posterior probability.

## Practical Tips

- Start with defaults if your BAM already carries `CB`, `UR`, `UY`, and `GX`.
- For separate samples or libraries that may reuse cell barcodes, use `--isolate-inputs` or include an existing source tag with `--source-tag`.
- Use `--no-gene` only if gene tags are absent or unreliable.
- Use `--end-bin` when you want molecule grouping to better respect transcript ends.
- Use `--emit-explain` when benchmarking or tuning parameters.
- Raise `--buckets` for larger datasets to spread memory across more buckets.
- Tighten `--min-merge-confidence` if you want more conservative UMI correction.

## Synthetic Truth Data And Tests

IsoUMI includes a deterministic synthetic truth generator for application-note
benchmarks and regression tests:

```bash
python3 scripts/generate_synthetic_truth.py --out-dir benchmarks/output/synthetic_truth
```

This writes:

- `synthetic_truth.sam`
- `synthetic_truth.truth.tsv`

The fixture covers UMI error correction within one junction context, same-UMI
different-junction separation, different-cell and different-chromosome
separation, non-spliced locus-bin correction, and missing-cell passthrough.

Run IsoUMI on the generated data:

```bash
src/isoumi \
  --bam benchmarks/output/synthetic_truth/synthetic_truth.sam \
  --out benchmarks/output/synthetic_truth/isoumi \
  --threads 1 \
  --buckets 8 \
  --ham 1 \
  --ratio 0.60 \
  --emit-tsv \
  --emit-explain
```

The full test suite now exercises both the hand-written smoke fixtures and the
synthetic truth generator:

```bash
make test
```

## Example Workflows

### Conservative correction

```bash
src/isoumi \
  --bam input.bam \
  --out conservative \
  --ham 1 \
  --ratio 0.05 \
  --min-merge-confidence 0.60 \
  --emit-explain
```

### Isoform-aware grouping

```bash
src/isoumi \
  --bam input.bam \
  --out isoform_run \
  --end-bin 25 \
  --emit-tsv \
  --emit-explain
```

### Tag remapping

```bash
src/isoumi \
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
- UMI quality is used for seed ranking and confidence scoring, not a full probabilistic error model.
- Final BAM linking depends on your local `htslib` setup; some static `htslib` builds may require extra libraries through `HTS_EXTRA_LIBS` or `CRYPTO_LIBS`.
- Output BAM is not coordinate-sorted by construction.

## Release Checks

Before tagging a publication release:

```bash
make clean
make
make test
make check-release
sh examples/run_minimal.sh
```

`make check-release` requires final publication metadata, including the public
repository URL and archive DOI, and intentionally fails while placeholders
remain.

The release checklist is maintained in [`docs/release_checklist.md`](./docs/release_checklist.md).

## Benchmark Scaffold

For publication benchmarks, use [`scripts/benchmark_isoumi.sh`](./scripts/benchmark_isoumi.sh) as the command template and document datasets under [`benchmarks/`](./benchmarks). The scaffold records dataset labels, input size and counts where available, exact commands, elapsed time, peak RSS, exit status, and parameters for default, transcript-end-aware, conservative correction, baseline, and optional sensitivity settings. Use [`scripts/summarize_benchmark_reports.py`](./scripts/summarize_benchmark_reports.py) to generate a manuscript-ready `summary.tsv`.

## Repository Layout

```text
.
├── CITATION.cff
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
├── Makefile
├── README.md
├── VERSION
├── benchmarks/
│   └── README.md
├── docs/
│   ├── application_note.md
│   └── release_checklist.md
├── examples/
│   ├── README.md
│   ├── minimal.sam
│   └── run_minimal.sh
├── tests/
│   └── smoke.sh
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

- `IsoUMI 0.1.0`
