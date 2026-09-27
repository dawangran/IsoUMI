# Contributing

Keep changes focused, reproducible, and documented. Report bugs with the
software version, command, error log, and a minimal input when possible.

## Development setup

Install a C11 compiler, GNU Make, HTSlib, pkg-config, Python 3, and samtools.
See the [installation instructions](README.md#installation) for package names.

```bash
make
make test
sh examples/run_minimal.sh
```

`make test` runs regression scripts covering UMI correction, synthetic
truth, alignment-direction grouping, BAM tags, header merging, non-primary
alignments, and input/output integrity, followed by release metadata checks.
Tests that inspect BAM records require samtools; a skipped test is not a
complete validation run. GitHub CI
installs all test dependencies and runs the full suite on Ubuntu.

## Changes and tests

- Document CLI behavior in `README.md` and detailed options in
  [`docs/parameters.md`](docs/parameters.md).
- Add or update regression tests when changing correction, grouping, report
  columns, or BAM semantics.
- Describe compatibility changes in `CHANGELOG.md`.
- Keep `VERSION`, `src/cli.h`, `CITATION.cff`, and documented versions consistent.
- Keep small, deterministic fixtures in `examples/` or the test scripts. Store
  generated analyses and manuscript materials outside version control, such as
  in the ignored `artifacts/` directory.

## Releases

Follow the [release checklist](docs/release_checklist.md), verify a clean build
and the full test suite, then tag the tested commit. Cite and archive the exact
software revision used for a published analysis.
