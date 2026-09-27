# Release checklist

Use a clean checkout of the commit intended for release.

## Version and metadata

- Keep `VERSION`, `src/cli.h`, `CITATION.cff`, `CHANGELOG.md`, and the README
  version in sync; confirm with `src/isoumi --version`.
- Confirm the repository URL, contributor metadata, and license.
- Update the changelog with user-visible changes and compatibility notes.
- Verify the documented container tag and WDL settings if publishing an image.

## Build and validation

Install all test dependencies, including Python 3 and samtools, then run:

```bash
make
make test
make check-release
sh examples/run_minimal.sh
```

Require successful GitHub CI for the release commit. Investigate any skipped
local tests before declaring validation complete. `make check-release` checks
version consistency, required metadata files, and unfinished repository or
archive placeholders; it does not verify remote availability or require a DOI.

## Publish and archive

- Create a version tag for the tested commit and describe its changes.
- Build container images from that commit and record their immutable digests.
- Archive the tagged source and reproducibility materials when preparing a
  published analysis; add a verified archive DOI to `CITATION.cff` when available.
- Record input accessions or checksums, reference/annotation versions, exact
  commands, and comparison settings with benchmark results.

Keep unpublished manuscripts and internal writing plans outside the software
repository. The [benchmark guide](../benchmarks/README.md) describes the supplied
fixtures and execution scripts.
