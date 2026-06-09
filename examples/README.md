# IsoUMI Minimal Example

This directory contains a tiny SAM fixture for checking command-line behavior and report semantics. It is not intended as a biological benchmark.

Build IsoUMI:

```bash
make
```

Run the example:

```bash
sh examples/run_minimal.sh
```

The script writes outputs under `examples/output/`:

- `minimal.dedup.bam`
- `minimal.molecules.tsv`
- `minimal.assignments.tsv`
- `minimal.corrections.tsv`

In this fixture, raw UMI `AAAT` is expected to collapse to `AAAA` in the first molecule group.
