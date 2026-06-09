# Contributing

IsoUMI is developed as a small command-line C tool. Keep changes focused, reproducible, and easy to audit.

## Development Setup

Build the command-line tool:

```bash
make
```

Run the smoke test:

```bash
make test
```

The smoke test constructs synthetic SAM inputs and verifies that IsoUMI produces a deduplicated BAM plus molecule, assignment, and correction reports.

## Code Guidelines

- Keep public command-line behavior documented in `README.md`.
- Add or update tests when changing UMI correction, grouping, report columns, or BAM tag semantics.
- Avoid changing output column names without documenting the compatibility impact.
- Keep the version in `VERSION`, `CITATION.cff`, `README.md`, and `src/cli.c` synchronized.

## Release Checklist

- `make clean && make`
- `make test`
- `src/isoumi --version`
- Confirm `README.md`, `CHANGELOG.md`, and `CITATION.cff` match the release version.
- Create a tagged release and archive it with a DOI provider such as Zenodo before manuscript submission.
