# IsoUMI Figure Design For Application Note

This document defines the figures for the IsoUMI Application Note. The main
figure should not be a generic workflow diagram. It should be an evidence-first
figure that shows why IsoUMI is needed and what it improves.

## Figure Strategy

Use one main figure in the article and place extended validation in the
supplement. Application Notes are short, so the main figure should combine the
problem, the method-specific advantage, and benchmark evidence.

Main figure:

- Figure 1: IsoUMI prevents transcript-context over-collapse while preserving
  explainable UMI correction.

Supplementary figures:

- Figure S1: Synthetic truth design and expected outcomes.
- Figure S2: Parameter sensitivity.
- Figure S3: Runtime and memory scaling.
- Figure S4: Example correction audit from `corrections.tsv`.

## Main Figure 1

### Title

IsoUMI preserves transcript-context-specific molecules during long-read
single-cell UMI deduplication.

### Layout

Use a 2 by 2 layout. Avoid a left-to-right workflow. Each panel should answer a
reviewer question.

```text
A. Why transcript context matters      B. Synthetic truth validation
C. Comparison with baselines           D. Real-data benchmark and auditability
```

## Figure 1A: The Over-Collapse Problem

### Goal

Show the biological/computational problem that motivates IsoUMI.

### Visual Design

Draw one cell and one gene with two long reads that have:

- same cell barcode
- same gene tag
- same raw UMI
- different splice-junction chains or different transcript ends

Represent the reads as exon blocks connected by intron arcs.

Use this contrast:

```text
Same CB + same GX + same UMI
Read 1: exon1--exon2--exon3
Read 2: exon1------exon3
```

On the right, show two outcomes:

```text
Gene-level UMI dedup: 1 molecule       over-collapse
IsoUMI: 2 molecule contexts            preserved
```

### Encoding

- Read 1: blue exon chain.
- Read 2: teal exon chain.
- Shared UMI: orange label.
- Incorrect merge: grey dashed bracket.
- Correct separation: two solid boxes.

### Panel Message

Reads that share cell, gene, and UMI can still represent different long-read
transcript contexts.

## Figure 1B: Synthetic Truth Validation

### Goal

Show that IsoUMI behaves correctly in controlled cases.

### Visual Design

Use a compact heatmap or checklist matrix.

Rows:

- UMI error within same junction
- Same UMI, different junction
- Same UMI, different cell
- Same UMI, different chromosome
- Non-spliced locus-bin correction
- Missing-cell passthrough

Columns:

- Expected
- IsoUMI default
- IsoUMI `--no-structure`
- UMI-tools

Cell values:

- check mark for expected behavior
- warning symbol for over-collapse or unsupported context
- dash for not applicable

### Example Matrix Logic

```text
Scenario                         IsoUMI default   --no-structure   UMI-tools
UMI error same junction          corrects         corrects         corrects
Same UMI different junction      separates        collapses        likely collapses*
Missing cell                     passthrough      passthrough      depends on setup
```

Use a footnote for UMI-tools:

```text
*UMI-tools is a general-purpose UMI deduplication baseline; exact behavior
depends on selected coordinate, cell, and gene options.
```

### Panel Message

The synthetic truth set isolates transcript-structure cases that are difficult
for non-structure-aware grouping.

## Figure 1C: Quantitative Comparison

### Goal

Compare IsoUMI against internal and external baselines.

### Visual Design

Use a grouped bar chart or dot plot.

X-axis:

- IsoUMI default
- IsoUMI `--end-bin 50`
- IsoUMI conservative
- IsoUMI `--no-structure`
- UMI-tools

Y-axis options:

- molecule count
- duplicate rate
- over-collapse events in synthetic truth

Best design: use two vertically stacked mini-plots in the same panel:

1. Molecule count retained.
2. Over-collapse count or discordant synthetic truth outcomes.

### Encoding

- IsoUMI modes: colored solid bars.
- Baselines: grey bars.
- UMI-tools: dark grey with outline.
- Mark the expected synthetic truth molecule count as a horizontal reference
  line.

### Panel Message

IsoUMI default preserves transcript-context-specific molecules better than
structure-agnostic baselines.

## Figure 1D: Real-Data Benchmark And Auditability

### Goal

Show practical usability on real or public long-read single-cell data and show
that correction decisions are inspectable.

### Visual Design

Use a split panel:

Left: benchmark dot/table inset.

```text
method | molecules | duplicate rate | corrected UMIs | runtime
```

Right: example correction table excerpt.

```text
raw_umi | corr_umi | count ratio | hamming | confidence | reason
AAAT    | AAAA     | 1/5         | 1       | 0.72       | count+quality
AATT    | AATT     | self        | 0       | 1.00       | self
```

### Encoding

- Runtime as small dots or a compact bar.
- Correction table as a clean inset with one highlighted accepted correction
  and one unchanged UMI.

### Panel Message

IsoUMI is practical to run and exposes correction decisions for manual review.

## Figure 1 Caption Draft

Figure 1. IsoUMI preserves transcript-context-specific molecules during
long-read single-cell UMI deduplication. (A) Reads from the same cell and gene
can share the same UMI while carrying distinct splice-junction chains or
transcript ends; structure-agnostic UMI deduplication can collapse these into a
single molecule, whereas IsoUMI keeps separate transcript contexts. (B)
Synthetic truth scenarios test within-context UMI correction, different-junction
separation, different-cell and different-chromosome separation, locus-bin
correction, and missing-cell passthrough. (C) IsoUMI modes are compared with a
no-structure internal baseline and UMI-tools, a general-purpose UMI
deduplication baseline. (D) Benchmark summaries and correction-explanation
records show runtime, molecule-level outputs, and inspectable UMI correction
decisions.

## Supplementary Figure S1: Synthetic Truth Design

Use one panel per synthetic scenario. This can be a clean table with small exon
sketches.

Required rows:

- `umi_error_same_junction`
- `same_umi_different_junction`
- `same_umi_different_cell`
- `same_umi_different_chromosome`
- `umi_error_locus_bin`
- `missing_cell_passthrough`

Columns:

- scenario
- read context
- raw UMI pattern
- expected IsoUMI result
- tested report

## Supplementary Figure S2: Parameter Sensitivity

Use line plots or dot plots.

Panels:

- molecule count vs `--ratio`
- corrected UMI count vs `--ratio`
- molecule count vs `--end-bin`
- correction count vs `--min-merge-confidence`

Recommended parameter grid:

```text
--ratio: 0.05, 0.10, 0.30, 0.60
--end-bin: off, 25, 50, 100
--min-merge-confidence: 0.00, 0.45, 0.60
```

## Supplementary Figure S3: Runtime And Memory

Use benchmark plots.

Panels:

- runtime vs read count or bucket count
- peak memory vs bucket count
- temporary disk usage vs bucket count

If peak memory is not available, report runtime and file sizes in a table
instead of plotting memory.

## Supplementary Figure S4: Correction Audit Examples

Show selected rows from `*.corrections.tsv`.

Include examples:

- accepted low-support UMI correction
- quality-supported correction
- unchanged self UMI
- correction blocked by direct Hamming threshold

## Data Needed To Finalize Figures

Minimum required:

- synthetic truth default IsoUMI output
- synthetic truth `--no-structure` output
- synthetic truth UMI-tools output
- one real/public long-read single-cell benchmark
- benchmark summary table from `scripts/summarize_benchmark_reports.py`

Optional but recommended:

- parameter sensitivity grid
- peak memory measurements
- one manually inspected high-depth cell/gene example
