# Changelog

## IsoUMI 0.1.0 - 2026-06-05

Initial publication-oriented release.

- Provides isoform-aware grouping from cell barcode, gene tag, strand, splice-junction or locus structure, and optional transcript-end bins.
- Corrects low-support UMIs toward higher-support UMIs using Hamming distance, count-ratio filtering, and optional quality-aware confidence scoring.
- Writes corrected UMI tags, duplicate flags, optional molecule tags, and optional TSV reports.
- Supports multi-BAM input with compatible reference dictionaries and metadata merging.
- Adds build/test entry points, a smoke test with synthetic SAM data, and publication metadata.
