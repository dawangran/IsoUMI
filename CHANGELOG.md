# Changelog

## Unreleased

- Allows otherwise matching gene-tagged reads expected to represent one molecule
  to remain together when they align in both directions. New
  `--strand-mode auto` (default) ignores BAM alignment
  direction for reads with a nonempty valid SAM `Z` string in the selected gene tag;
  missing/unusable gene tags and `--no-gene` retain direction separation.
- Adds `--strand-mode alignment` to reproduce historical grouping and
  `--strand-mode ignore` to omit direction regardless of gene annotations.
  The policy applies to exact, compatible, and no-structure grouping with both
  correction methods. Ignored direction uses genomic `EL`/`ER` end bins.
  Default grouping keys and molecule IDs change for gene-tagged reads even
  when molecule counts do not. Raw UMI strings, sequences, and alignment orientation
  remain unchanged; transcript strand is not inferred.
  Mixed directions depend on library preparation and upstream orientation;
  `auto` does not detect either, and the regression evidence is synthetic.
- Adds experimental `--structure-mode compatible` with direct splice-coordinate
  tolerance and truncation-aware structural assignment. Raw UMI Hamming
  neighborhoods limit the search; UMI correction is recomputed within final
  structure groups. Exact grid grouping remains the default.
- Retains conflicting low-support structures and reports ambiguous structural
  assignments instead of allowing partial reads to join incompatible anchors.

- Adds `--strip-pg` to consistently omit header `@PG` records, record-level
  `PG` tags, and IsoUMI's own program record from temporary and final BAMs.
- Stream-aggregates repeated grouping-key/UMI observations, avoiding per-read
  strings and full-read sorting in each bucket.
- Reuses the aggregation pass for molecule representative selection and
  non-primary discovery, reducing each bucket from four or five full scans to
  two scans without non-primary records or three scans with them.
- Uses fast compression for transient bucket BAMs and applies `--threads` to
  input decoding and final BAM compression.
- Reclaims source buckets after successful per-bucket deduplication, then
  deduplicated BAM and report buckets after successful final-output flushing;
  `--keep-tmp` preserves all intermediate generations for debugging.

## IsoUMI 0.1.1 - 2026-08-27

- Preserves multi-input program provenance by namespacing later `@PG` chains
  and rewriting header `ID`/`PP` plus record-level `PG:Z:` references together;
  output headers also record the IsoUMI version and command line.
- Counts only primary alignments for UMI correction and molecule selection;
  secondary and supplementary records inherit the matching primary read's
  corrected UMI and duplicate status.
- Rejects conflicting reference `M5`/`AS` metadata and conflicting or invalid
  SAM tag configurations before processing.
- Adds optional UMI-tools-style directional UMI correction with the `2n-1` count rule and transitive graph grouping.
- Extends correction reports with method, immediate-parent, edge-distance, and path-length audit fields.
- Preserves direct count-ratio correction as the backward-compatible default.

## IsoUMI 0.1.0 - 2026-06-05

Initial publication-oriented release.

- Provides isoform-aware grouping from cell barcode, gene tag, strand, splice-junction or locus structure, and optional transcript-end bins.
- Corrects low-support UMIs toward higher-support UMIs using Hamming distance, count-ratio filtering, and optional quality-aware confidence scoring.
- Writes corrected UMI tags, duplicate flags, optional molecule tags, and optional TSV reports.
- Supports multi-BAM input with compatible reference dictionaries and metadata merging.
- Adds build/test entry points, a smoke test with synthetic SAM data, and publication metadata.
