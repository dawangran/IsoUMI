# IsoUMI Benchmark Scaffold

This directory documents the benchmark runs expected for a publication release. It intentionally does not include large input data.

Use `scripts/benchmark_isoumi.sh` with a real or simulated long-read single-cell BAM:

```bash
sh scripts/benchmark_isoumi.sh path/to/input.bam benchmarks/output
```

For the synthetic truth benchmark, generate a deterministic SAM fixture and truth table:

```bash
python3 scripts/generate_synthetic_truth.py --out-dir benchmarks/output/synthetic_truth
src/isoumi \
  --bam benchmarks/output/synthetic_truth/synthetic_truth.sam \
  --out benchmarks/output/synthetic_truth/isoumi \
  --tmp-dir benchmarks/output/synthetic_truth/tmp \
  --threads 1 \
  --buckets 8 \
  --ham 1 \
  --ratio 0.60 \
  --emit-tsv \
  --emit-explain
```

The generated truth table marks the scenarios needed for the Application Note: within-junction UMI error correction, same-UMI different-junction over-collapse under a naive baseline, different-cell/chromosome separation, locus-bin correction, and missing-cell passthrough.

The script runs three representative configurations and writes a compact `metrics.tsv` file:

- `default`: default UMI correction settings.
- `endbin50`: transcript-end-aware grouping with `--end-bin 50`.
- `conservative`: stricter count-ratio and merge-confidence thresholds.

For a manuscript, report:

- input dataset name and accession or simulation parameters;
- number of reads, cells, and tagged UMIs;
- command line and IsoUMI version;
- runtime and peak memory where available;
- output molecule counts and duplicate rates;
- comparison baseline and its exact command line.
