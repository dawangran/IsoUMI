# IsoUMI Application Note Manuscript Plan

This document is a working plan for preparing an Application Note submission.
It assumes a Bioinformatics-style short software article: approximately 2600
words, or about 2000 words plus one main figure.

## Core Message

IsoUMI provides isoform-aware UMI correction and molecule-level deduplication
for long-read single-cell transcriptomic BAM files. Its main contribution is
that UMI correction is constrained within cell and transcript context rather
than only cell/gene context, reducing over-collapse of reads from distinct
splice-junction or transcript-end structures.

## Manuscript Structure

### Title

IsoUMI: isoform-aware UMI correction and molecule-level deduplication for
long-read single-cell transcriptomics

### Abstract

Target length: 150-200 words.

Required content:

- Problem: long-read single-cell isoform analysis requires molecule-level
  deduplication that respects transcript structure.
- Method: IsoUMI groups reads by cell barcode, gene, strand, splice-junction or
  locus context, and optional transcript-end bins; it corrects low-support UMIs
  using direct seed-based sequence, count-ratio, and optional quality filters.
- Output: corrected BAM plus molecule, assignment, and correction-explanation
  tables.
- Availability: GitHub URL, archived DOI, license, test data.

### 1. Introduction

Target length: 350-450 words.

Key points:

- UMIs are essential for molecule-level deduplication in single-cell assays.
- Long-read single-cell data adds isoform-resolved information, but naive UMI
  deduplication can over-collapse reads from distinct transcript structures.
- Existing short-read assumptions are not always appropriate for long-read
  transcriptomic alignments because splice-junction chains and transcript ends
  matter.
- IsoUMI addresses this by making transcript context part of the UMI correction
  unit.

### 2. Implementation

Target length: 550-700 words.

Describe:

- Input expectations: SAM/BAM, cell barcode, raw UMI, optional UMI quality, gene
  tags.
- Cell-bucket sharding and multi-input header compatibility checks.
- Grouping key:
  - cell barcode
  - gene or `NA`
  - reference target
  - strand
  - splice-junction hash for spliced reads
  - locus-bin hash for non-spliced reads
  - optional 5'/3' transcript-end bins
- `--no-structure` baseline for evaluating the contribution of transcript
  context.
- Seed-based UMI correction:
  - strongest seed UMIs processed first
  - direct raw-to-seed Hamming, count-ratio, and optional confidence filters
  - no chain-based correction beyond the specified direct distance
- Representative read selection:
  - mapping quality
  - aligned reference span
  - query length
  - stable input order
- Output tags and reports.

### 3. Results

Target length: 650-850 words.

Recommended subsections:

1. Synthetic truth validation
   - Use `scripts/generate_synthetic_truth.py`.
   - Report expected corrections and pass-through behavior.
   - Show default IsoUMI keeps same-UMI different-junction molecules separate.
   - Show `--no-structure` baseline over-collapses the same synthetic case.

2. Public long-read single-cell dataset benchmark
   - Report dataset source and accession.
   - Run `scripts/benchmark_isoumi.sh`.
   - Compare default, transcript-end-aware, conservative, and `--no-structure`
     modes.
   - Report runtime, memory if available, molecule counts, duplicate rates, and
     correction counts.

3. Parameter sensitivity
   - Summarize molecule count stability across `--ham`, `--ratio`, and
     `--end-bin`.
   - Use `--emit-explain` to inspect corrections with low confidence.

### 4. Availability And Implementation

Target length: 100-150 words.

Include:

- Source code URL.
- Archived release DOI.
- Version.
- License.
- Platform and dependencies.
- Test command: `make test`.
- Synthetic truth generation command.

### 5. Limitations

Target length: 120-180 words.

State clearly:

- UMI sequence comparison currently uses Hamming distance.
- Merge confidence is heuristic, not a posterior probability.
- Transcript-end bins are useful but depend on alignment and truncation quality.
- Output BAM is bucket-concatenated and should be treated as unsorted unless
  sorted downstream.

## Main Figure

Use one evidence-first main figure, not a generic workflow diagram. The figure
should show the over-collapse problem, the synthetic truth validation, the
baseline comparison, and real-data usability. Detailed design notes are in
`docs/figure_design.md`.

### Figure 1A: Why Transcript Context Matters

Show two reads from the same cell and gene with the same UMI but different
splice-junction chains or transcript ends.

Message: structure-agnostic UMI deduplication can over-collapse distinct
long-read transcript molecules.

### Figure 1B: Synthetic Truth Validation

Show a matrix of expected outcomes across synthetic truth scenarios for IsoUMI
default, IsoUMI `--no-structure`, and UMI-tools.

Message: synthetic data isolate the cases where transcript context matters.

### Figure 1C: Quantitative Baseline Comparison

Show molecule counts, duplicate rates, and over-collapse or discordant outcomes
for IsoUMI modes, `--no-structure`, and UMI-tools.

Message: IsoUMI preserves transcript-context-specific molecules better than
structure-agnostic baselines.

### Figure 1D: Real-Data Benchmark And Auditability

Show a compact real-data benchmark table or dot plot together with selected
`corrections.tsv` rows.

Message: IsoUMI is practical to run and its correction decisions are
inspectable.

## Main Table

Table 1 should summarize benchmark results.

Suggested columns:

- dataset
- mode
- reads
- cells
- molecules
- duplicate rate
- corrected UMI count
- runtime
- peak memory

Rows:

- synthetic truth, default
- synthetic truth, `--no-structure`
- real dataset, default
- real dataset, `--end-bin 50`
- real dataset, conservative
- real dataset, `--no-structure`

## Supplementary Materials

### Table S1: Parameter Guide

Use `docs/parameters.md` as the source. Include parameter name, default, stage,
and effect.

### Table S2: Synthetic Truth Scenarios

Summarize each scenario generated by `scripts/generate_synthetic_truth.py`.

Suggested columns:

- scenario
- expected behavior
- tested output file
- tested condition

### Figure S1: Parameter Sensitivity

Plot molecule counts and correction counts across:

- `--ratio 0.05`, `0.10`, `0.30`, `0.60`
- `--end-bin 25`, `50`, `100`
- `--min-merge-confidence 0.00`, `0.45`, `0.60`

### Figure S2: Runtime And Memory Scaling

Plot runtime and peak memory across read-count subsets or bucket counts.

## Required Before Submission

- Replace placeholder GitHub URL in `README.md`, `CITATION.cff`, and manuscript.
- Tag a release and archive it in Zenodo or another stable repository.
- Add the archived DOI to the manuscript availability section.
- Run `make test` on a clean machine or container.
- Run `scripts/benchmark_isoumi.sh` on at least one public or shareable real
  long-read single-cell dataset.
- Generate `summary.tsv` with `scripts/summarize_benchmark_reports.py`.
- Prepare Figure 1 as vector artwork.
- Prepare one supplementary PDF with parameter table, synthetic truth design,
  benchmark commands, and extended results.
- Add an AI assistance disclosure if required by the target journal.
