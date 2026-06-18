# IsoUMI Release Checklist

Use this checklist before tagging a release for manuscript submission.

## Required Metadata

- `VERSION` contains the release version.
- `src/cli.c` prints the same version through `src/isoumi --version`.
- `README.md`, `CHANGELOG.md`, and `CITATION.cff` use the same release version.
- `CITATION.cff` contains the final repository URL and archive DOI.
- `LICENSE` matches the approved project license.

## Build And Test

```bash
make clean
make
make test
make check-release
sh examples/run_minimal.sh
```

`make check-release` runs the strict metadata gate and fails until placeholder
repository URLs and archive DOI text have been replaced.

## Manuscript Artifacts

- Application note draft updated in `docs/application_note.md`.
- Workflow figure exported and referenced in the manuscript.
- Benchmark scripts and datasets archived or linked; use `scripts/benchmark_isoumi.sh` as the command template.
- Release tag created.
- DOI archived through Zenodo or another repository.

## Submission Notes

For an Application Note submission, report the exact release version, repository URL, license, dependencies, and availability of example or benchmark data.
