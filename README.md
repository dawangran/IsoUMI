<p align="center"><img src="docs/assets/isoumi-icon.png" width="112" alt="IsoUMI icon" /></p>
<h1 align="center">IsoUMI</h1>
<p align="center">Alignment-context-aware UMI correction for long-read single-cell sequencing</p>
<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue" alt="License: MIT" /></a>
  <a href="src/Makefile"><img src="https://img.shields.io/badge/language-C11-555555" alt="Language: C11" /></a>
  <a href="VERSION"><img src="https://img.shields.io/badge/version-0.1.1-168C8C" alt="Source version: 0.1.1" /></a>
  <a href="https://github.com/dawangran/IsoUMI/actions/workflows/ci.yml"><img src="https://github.com/dawangran/IsoUMI/actions/workflows/ci.yml/badge.svg" alt="Build and tests" /></a>
</p>
<p align="center">
  <a href="#quick-start">Quick start</a> ·
  <a href="docs/parameters.md">Parameter guide</a> ·
  <a href="wdl/README.md">WDL workflow</a> ·
  <a href="https://github.com/dawangran/IsoUMI/issues">Issues</a>
</p>

IsoUMI corrects unique molecular identifiers (UMIs) and marks duplicate reads in
long-read single-cell BAM files. It restricts UMI comparisons by cell barcode,
gene, configurable alignment direction, and alignment context: splice-junction
chains for spliced reads, or genomic start/end bins for non-spliced reads.
Optional alignment-end bins and sample boundaries provide additional control
over grouping.

The tool writes corrected UMI tags, duplicate labels, and optional audit tables.
It operates on existing alignments and annotations; it does not align reads,
assign genes, or reconstruct transcript isoforms.

## Installation

The current source version is **0.1.1**. Changes under
[Unreleased](CHANGELOG.md#unreleased), including the strand-grouping fix, require
building the current source and are not part of the published 0.1.1 release.
Build from source on Linux or macOS, or use the
[Docker image](#docker-and-jupyterlab) for the published release.

### Build from source

Requirements: a C11 compiler, GNU Make, HTSlib and its compression/network
libraries. `pkg-config` is recommended for locating HTSlib. Install Python 3
and samtools for the full test suite; samtools also provides BAM filtering and
indexing.

On Ubuntu/Debian, install build and test dependencies:

```bash
sudo apt-get update
sudo apt-get install --yes --no-install-recommends \
  build-essential pkg-config libhts-dev samtools python3
```

The package manager installs the compression and network libraries required by
`libhts-dev`; a separate curl development package is not needed.

Then build and check the executable:

```bash
git clone https://github.com/dawangran/IsoUMI.git
cd IsoUMI
make
src/isoumi --version
make test
```

The executable is `src/isoumi`. The default build avoids CPU-specific
instructions; `make NATIVE=1` enables optimization for the build machine.

For an HTSlib installation outside standard paths:

```bash
make -C src HTSLIB_INC=/path/to/htslib/include HTSLIB_LIB=/path/to/htslib/lib
```

Additional linking requirements can be supplied through `HTS_EXTRA_LIBS` or
`CRYPTO_LIBS`, for example `make CRYPTO_LIBS="-lcrypto"`.

**Threading:** Linux builds enable OpenMP by default. On macOS, OpenMP is off by
default; `--threads` still controls HTSlib input decoding and final BAM
compression, but bucket processing is serial and prints a runtime notice.
Configure `OPENMP_CFLAGS` and `OPENMP_LDFLAGS` for an OpenMP-capable toolchain to
enable parallel bucket processing on macOS.

## Quick start

Run on a BAM containing `CB` cell barcodes and `UR` raw UMIs:

```bash
src/isoumi \
  --bam input.bam \
  --out sample \
  --emit-tsv \
  --emit-explain
```

This uses **ratio correction** with Hamming distance `1` and a maximum
low-count/high-count ratio of `0.10`. In current source, the default
`--strand-mode auto` lets opposite alignment directions share a molecule group
when the selected gene tag is present and nonempty. It produces:

| File | Contents |
| --- | --- |
| `sample.dedup.bam` | All input records with corrected UMI and duplicate annotations |
| `sample.molecules.tsv` | Corrected UMI groups and their primary-read support |
| `sample.assignments.tsv` | Eligible primary-read assignments and duplicate labels |
| `sample.corrections.tsv` | Raw-to-corrected UMI mappings and correction diagnostics |

Despite its `.dedup.bam` suffix, the BAM **retains duplicate records**. To remove
records marked as duplicates, then coordinate-sort and index the result:

```bash
samtools view -b -F 1024 sample.dedup.bam \
  | samtools sort -o sample.filtered.sorted.bam
samtools index sample.filtered.sorted.bam
```

This filter only removes records carrying the SAM duplicate bit. Unmapped,
secondary, supplementary, and other passthrough records can remain; it is not a
strict one-record-per-molecule filter. IsoUMI's BAM output itself is unsorted
(`SO:unknown`).

For a small, self-contained example:

```bash
sh examples/run_minimal.sh
```

The [example](examples/README.md) writes to `examples/output/` and demonstrates
correction of `AAAT` to `AAAA` in one UMI group.

## Input requirements

IsoUMI accepts one or more SAM/BAM files without requiring pre-sorting or an
index. Preprocessing must supply the cell and UMI tags used by your analysis;
tag names are configurable.

| Tag | Default | Role |
| --- | --- | --- |
| Cell barcode | `CB` | Required for a mapped primary read to participate in correction |
| Raw UMI | `UR` | Required for UMI correction and duplicate grouping |
| Gene assignment | `GX` | Included in grouping by default; missing values form an `NA` group |
| UMI quality | `UY` | Optional Phred+33 string, with the same length as the UMI |

Use `--no-gene` to omit gene annotations, or remap tags such as
`--cell-tag XC --umi-tag XM --gene-tag GN`. Tag names must be valid, distinct
SAM tags; `PG` is reserved for program provenance. UMI comparison uses Hamming
distance and therefore compares equal-length strings.

The default `--strand-mode auto` trusts a nonempty SAM `Z` string in the selected gene
tag and ignores the read's BAM alignment direction for that record. Reads with
the same cell, gene, UMI, and alignment context can therefore form one molecule
even when their alignments use both directions. Reads without a usable gene tag,
or all reads with `--no-gene`, retain alignment-direction separation. Use
`--strand-mode alignment` to reproduce the earlier direction-based grouping, or
`--strand-mode ignore` to ignore direction even without gene annotations.

Whether opposite alignment directions can represent one RNA molecule depends
on library preparation and upstream read orientation. `auto` checks gene tags;
it does not detect the protocol or whether reads have already been oriented.
Choose `alignment` when direction remains an informative molecule boundary.

This changes grouping keys and molecule IDs for gene-tagged reads, even when
molecule counts do not change. IsoUMI does not infer RNA transcription strand,
reverse-complement UMI tags, or change alignment orientation. Gene annotations
must distinguish features that should remain separate: same-gene antisense
reads cannot be resolved by this policy alone, and UMI collisions remain possible.
See the [strand-mode guide](docs/parameters.md#--strand-mode-autoalignmentignore-default-auto).

### Multiple inputs

Technical splits from the same barcode namespace can be processed together:

```bash
src/isoumi --bam lane1.bam --bam lane2.bam --out combined --emit-tsv
```

Alternatively, use `--bam-list bam_files.txt`, with one path per line. Blank
lines and lines beginning with `#` are ignored. Inputs are deduplicated together
by default. Use `--isolate-inputs` for independent samples that may reuse cell
barcodes, or `--source-tag RG` when read-group tags define the library boundary.

Reference dictionaries must agree in target name, order, and length.
Conflicting reference metadata or read-group definitions are rejected. Program
headers and record-level `PG` tags are merged consistently, and IsoUMI records
its version and command line in a new `@PG` entry. See the
[parameter guide](docs/parameters.md#inputs-and-outputs) for header validation
and the optional `--strip-pg` behavior.

## Output semantics

- **Corrected UMI:** written to `UB` by default; change it with `--umi-out`.
- **Duplicate status:** `DA:i:1` marks a duplicate; `DA:i:0` marks a
  representative or passthrough record. The SAM duplicate bit (`0x400`) is set
  or cleared consistently unless `--no-bam-dup-flag` is used.
- **Representative selection:** higher mapping quality, then longer aligned
  reference span, then longer query length, then earlier input order.
- **Molecule IDs:** `--mol-tag` writes the group and corrected UMI identifier to
  a BAM tag. IDs also appear in the assignment report.

Only eligible primary alignments contribute to UMI counts and molecule reports.
Unmapped reads and mapped primary reads without a cell barcode or raw UMI pass
through with duplicate status `0`. Their corrected UMI equals the raw UMI when
available; if the raw UMI is absent, stale corrected-UMI and molecule-ID tags are
removed.

Secondary and supplementary alignments inherit corrected UMI and duplicate
status from an unambiguous matching primary record. They do not inflate UMI
counts or assignment reports. Mapped non-primary records must carry the same
cell barcode as their primary; those without a cell tag are rejected. Records
without a matching eligible primary pass through without correction. Ambiguous
primary read identities are rejected; use `--isolate-inputs` when independent
files may reuse read names.

## Method

IsoUMI partitions records into cell-barcode buckets, then corrects UMIs within
keys containing cell, gene (unless disabled), optional source, reference,
alignment direction according to `--strand-mode`, and alignment context.

| Read context | Default grouping |
| --- | --- |
| Spliced (`N` in CIGAR) | Junction chain with boundaries rounded to the nearest 10 bp |
| Non-spliced | Pair of genomic start and exclusive end bins, each 1,000 bp wide |
| Optional alignment ends | Genomic left/right bins when direction is ignored; direction-oriented 5′/3′ bins otherwise, set with `--end-bin` |

`--sj-jitter` specifies rounding to a fixed grid, not a pairwise distance
threshold. Reads on opposite sides of a rounding boundary can fall into
different groups even when their coordinates differ by only one base.
`--no-structure` omits junction, locus, and end terms while retaining cell,
gene, source, reference, and the selected alignment-direction policy. The same
`--strand-mode` policy applies to exact, compatible, and no-structure grouping,
with either correction method.

Two correction methods are available:

- **`ratio` (default):** each raw UMI must satisfy the Hamming-distance,
  count-ratio, and optional confidence thresholds directly against its final
  seed. Seed priority is count, optional mean UMI quality, then lexicographic
  order.
- **`directional`:** an edge from UMI A to B requires
  `count(A) >= 2 × count(B) − 1` and distance within `--ham`. Reachable UMIs are
  assigned to a seed in descending count and lexicographic order. Transitive
  paths can yield a final raw-to-seed distance greater than `--ham`.

Directional mode uses the UMI-tools count/network rule within IsoUMI's grouping
contexts. Different grouping definitions can produce different results from
standalone UMI-tools. In directional mode, `--ratio` is ignored, quality does
not determine network membership, and positive `--min-merge-confidence` values
are rejected.

**Interpretation:** alignment context can separate UMIs associated with
different transcript structures, but alignment variation and incomplete reads
can also split evidence from the same molecule. More reported groups do not,
on their own, establish greater molecular accuracy. Reported merge confidence
is a heuristic score, not a posterior probability. Directional pairwise
comparison can be costly for groups with many distinct UMIs.

## Compatible structural grouping (experimental)

The optional compatible mode groups reads affected by truncation or small
splice-coordinate shifts:

```bash
src/isoumi --bam input.bam --out compatible \
  --structure-mode compatible --sj-tolerance 10 \
  --min-structure-support 3 \
  --correction-method directional --ham 1 --emit-tsv --emit-explain
```

This mode compares observed alignment geometry directly. It constructs local
UMI candidate neighborhoods, groups compatible structures, and recomputes UMI
correction from the raw counts within each final structure group. Partial reads
can attach to compatible longer structures without joining conflicting anchors.
Low-support observations with incompatible geometry remain separate.

`--structure-mode exact` preserves the existing grid-based grouping and is the
default. The compatible mode is incompatible with `--no-structure` and
`--end-bin`; `--sj-jitter` and `--locus-bin` configure only exact mode. See the
[parameter guide](docs/parameters.md#--structure-mode-exactcompatible-default-exact).

Compatibility is an inference from alignments, not molecular ground truth.
Shared truncated fragments can be ambiguous, nearby genuine splice sites can
fall within the tolerance, and overlapping unspliced reads may not distinguish
alternative transcript ends. Evaluate assignment ambiguity, splitting and
mixing alongside molecule counts before using this mode for quantification.

## Common configurations

These examples illustrate specific analysis choices; they are not universal
parameter recommendations. See the [parameter guide](docs/parameters.md) for
all options, defaults, and interactions.

```bash
# Directional correction within alignment-context groups
src/isoumi --bam input.bam --out directional \
  --correction-method directional --ham 1 --emit-tsv --emit-explain

# Stricter direct ratio correction
src/isoumi --bam input.bam --out conservative \
  --ratio 0.05 --min-merge-confidence 0.60 --emit-explain

# Add 50 bp alignment-end bins
src/isoumi --bam input.bam --out end_aware --end-bin 50 --emit-tsv

# Reproduce historical alignment-direction grouping
src/isoumi --bam input.bam --out legacy \
  --strand-mode alignment --emit-tsv

# Compare the same correction method without alignment structure
src/isoumi --bam input.bam --out no_structure \
  --correction-method directional --no-structure --emit-tsv --emit-explain

# Keep independent inputs separate even when barcodes are reused
src/isoumi --bam sample1.bam --bam sample2.bam --out separate \
  --isolate-inputs --emit-tsv
```

For large inputs, `--threads`, `--buckets`, and `--tmp-dir` control parallelism
and scratch storage. Use a new temporary directory on local SSD/NVMe storage
and allow space for the growing final output. More buckets can reduce memory
per bucket, although one high-complexity cell still occupies a single bucket.
Temporary files are reclaimed incrementally unless `--keep-tmp` is set.

## Docker and JupyterLab

The container includes IsoUMI, samtools, Python, and JupyterLab:

```bash
docker pull dawang02/isoumi:0.1.1
docker run --rm -v "$PWD:/data" dawang02/isoumi:0.1.1 \
  isoumi --bam /data/input.bam --out /data/sample --threads 8
```

Start JupyterLab with a private token and persist work in the current directory:

```bash
docker run --rm -p 127.0.0.1:8888:8888 \
  -e JUPYTER_TOKEN='replace-with-a-private-token' \
  -v "$PWD:/workspace" dawang02/isoumi:0.1.1
```

Open `http://localhost:8888` with the configured token. If the token variable is
omitted, Jupyter generates a token and prints the access URL in the container
logs. The [quick-start notebook](notebooks/IsoUMI_quickstart.ipynb) provides an
interactive example.

The [WDL workflow](wdl/README.md) runs directional correction without gene
grouping or quality-aware scoring and produces both the annotated BAM and a
filtered, sorted, indexed BAM. Its analysis settings differ from CLI defaults;
review them before running it with miniwdl or Cromwell.

## Reproducibility and citation

Run `make test` for regression checks and follow the
[benchmark guide](benchmarks/README.md) for synthetic fixtures, benchmark
commands, and report summaries. Small fixtures test defined behaviors; they do
not establish performance across biological datasets.

For reproducible analyses, record the source commit and version, container
digest if applicable, exact command, input accession or checksum, reference
and annotation versions, and baseline settings. Benchmark results apply to the
revision and configuration tested; rerun them when either changes.

Software citation metadata is available in [`CITATION.cff`](CITATION.cff).
Cite the exact version or commit used. See the
[release checklist](docs/release_checklist.md) for publication and archival
requirements.

IsoUMI is distributed under the [MIT license](LICENSE). For bug reports or
feature requests, open a [GitHub issue](https://github.com/dawangran/IsoUMI/issues)
with the version, command, error log, and a minimal example when possible.
