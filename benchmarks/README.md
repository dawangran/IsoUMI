# IsoUMI Benchmark Scaffold

This directory documents the benchmark runs expected for a publication release. It intentionally does not include large input data.

Use `scripts/benchmark_isoumi.sh` with a real or simulated long-read single-cell BAM:

```bash
sh scripts/benchmark_isoumi.sh \
  --dataset <dataset_name_or_accession> \
  --threads 8 \
  --buckets 128 \
  path/to/input.bam \
  benchmarks/output
```

For the synthetic truth benchmark, generate a deterministic SAM fixture and truth table:

```bash
python3 scripts/generate_synthetic_truth.py --out-dir benchmarks/output/synthetic_truth
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

The generated truth table marks the scenarios needed for the Application Note: within-junction UMI error correction, same-UMI different-junction over-collapse under a naive baseline, different-cell/chromosome separation, locus-bin correction, and missing-cell passthrough.

The same generator is exercised by the regression suite:

```bash
sh tests/synthetic_truth.sh src/isoumi
```

The script runs representative configurations and writes `metrics.tsv`,
per-mode logs, and IsoUMI report files:

- `default`: default UMI correction settings.
- `endbin50`: transcript-end-aware grouping with `--end-bin 50`.
- `conservative`: stricter count-ratio and merge-confidence thresholds.
- `baseline_no_structure`: baseline grouping with `--no-structure`, omitting splice-junction, locus, and transcript-end context.

`metrics.tsv` records the fields needed for manuscript tables:

- dataset label, IsoUMI version, input path, input format, and input size;
- input read, cell-barcode, and UMI counts when available;
- threads, buckets, start/end UTC timestamps, elapsed seconds, and peak RSS;
- exit status, mode-specific parameters, and exact command line.

For parameter-sensitivity runs, use:

```bash
sh scripts/benchmark_isoumi.sh \
  --dataset <dataset_name_or_accession> \
  --case-set sensitivity \
  path/to/input.bam \
  benchmarks/output_sensitivity
```

Use `--case-set all` to run both the representative benchmark and sensitivity
grid in one output directory.

After the benchmark finishes, summarize report-level metrics for manuscript
tables:

```bash
python3 scripts/summarize_benchmark_reports.py benchmarks/output
```

This writes `benchmarks/output/summary.tsv` with input metadata, assigned read
and cell counts, molecule counts, molecule-support totals, duplicate rates,
correction counts, correction rates, mean merge confidence, runtime, peak RSS,
version, parameters, and exact command line. Missing report files are treated as
errors by default; use `--allow-missing` only for exploratory debugging.

For a manuscript, report:

- input dataset name and accession or simulation parameters;
- number of reads, cells, and tagged UMIs;
- command line and IsoUMI version;
- runtime and peak memory where available;
- output molecule counts and duplicate rates;
- comparison baseline and its exact command line.
