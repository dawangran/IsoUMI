# IsoUMI Application Note Draft

## Title

IsoUMI: isoform-aware UMI correction and molecule-level deduplication for long-read single-cell transcriptomics

## Abstract

Long-read single-cell transcriptomic assays enable isoform-resolved analysis but require molecule-level deduplication strategies that respect transcript structure. IsoUMI is a command-line tool that corrects raw UMI tags and marks duplicates in long-read single-cell BAM files. It groups reads by cell barcode, optional gene tag, strand, splice-junction or locus structure, and optional transcript-end bins, then collapses low-support UMIs toward higher-support UMIs using sequence distance, count-ratio filtering, and optional UMI quality information. IsoUMI outputs a corrected BAM and optional molecule, assignment, and correction-explanation reports, supporting reproducible inspection of UMI correction decisions.

## Availability

IsoUMI is implemented in C and uses htslib for BAM/SAM I/O. The software is released under the MIT license.

- Source code: https://github.com/your-org/IsoUMI
- Version described here: 0.1.0
- Archive DOI: to be added after release

## Implementation

IsoUMI runs in three main phases. First, input reads are split into cell-barcode buckets so independent cells can be processed in parallel. Second, each bucket is scanned to build grouping keys from cell, optional source or input scope, and transcript context. For spliced reads, the key includes a hash of splice-junction coordinates after optional boundary jitter rounding; for non-spliced reads, the key includes genomic locus bins. Optional transcript-end bins can further separate molecule contexts. Third, unique UMIs are counted within each group, candidate seed UMIs are processed from strongest to weakest, and lower-support UMIs are corrected only when they directly satisfy the Hamming-distance, count-ratio, and optional confidence filters against that seed.

Corrected UMI values are written to the output BAM tag selected by `--umi-out`; duplicates are marked with both the tag selected by `--dup-flag` and, by default, the standard SAM duplicate bit. Within each corrected molecule, the retained representative read is selected by mapping quality, aligned reference span, query length, and stable input order. Unmapped reads and mapped reads without a cell barcode are passed through and do not participate in correction. Optional reports summarize corrected molecules, read assignments, and correction decisions.

## Use Case

A typical run starts from one or more long-read single-cell BAM files containing cell barcode (`CB`) and raw UMI (`UR`) tags:

```bash
src/isoumi \
  --bam input.bam \
  --out sample \
  --threads 8 \
  --buckets 128 \
  --emit-tsv \
  --emit-explain \
  --end-bin 50
```

## Evaluation Plan

The publication version should include the following evaluations before submission:

- Runtime and memory usage across increasing read counts and bucket settings.
- Molecule count stability under different `--ham`, `--ratio`, and `--end-bin` values.
- Comparison with `--no-structure`, a baseline UMI deduplication mode that does not use splice-junction, locus, or transcript-end context.
- Manual inspection of correction reports on representative high-depth cells.
- Reproducible scripts and a small public or simulated dataset. The repository includes `scripts/generate_synthetic_truth.py` and `tests/synthetic_truth.sh` for deterministic synthetic-data generation and regression testing.

## Figure Plan

Figure 1 should show the IsoUMI workflow:

1. Tagged long-read single-cell BAM input.
2. Cell barcode sharding.
3. Transcript-context grouping from gene, strand, splice-junction/locus hash, and optional transcript-end bins.
4. Seed-based UMI correction using count ratio, Hamming distance, and optional quality-aware confidence.
5. Representative read selection within corrected molecules.
6. Corrected BAM and TSV reports.

## Current Limitations

IsoUMI currently uses Hamming distance for UMI sequence comparison and does not implement a full probabilistic sequencing-error model. Merge confidence values are heuristic scores rather than posterior probabilities. Output BAM records are bucket-concatenated and should be treated as unsorted unless sorted downstream.
