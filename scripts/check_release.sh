#!/bin/sh
set -eu

version=$(cat VERSION)
expected_name="IsoUMI"
expected_version="$expected_name $version"

grep -Fq "#define ISOUMI_VERSION_NUMBER \"$version\"" src/cli.h
grep -Fq "The current source version is **$version**." README.md
grep -Fq "version-$version-" README.md
grep -q "version: \"$version\"" CITATION.cff
grep -q "## $expected_name $version" CHANGELOG.md
grep -q "BIN = isoumi" src/Makefile
grep -q "src/isoumi" .gitignore
grep -q -- '--correction-method' README.md
grep -Fq 'seed_count >= 2 * raw_count - 1' docs/parameters.md
grep -q 'parent_umi' docs/parameters.md

if [ -x src/isoumi ]; then
  actual=$(src/isoumi --version 2>&1)
  test "$actual" = "$expected_version"
fi

stale_name_upper="$(printf '%s%s' Long UMI)"
stale_name_lower="$(printf '%s%s' long umi)"
stale_version_prefix='0[.]3[.]'
text_paths='README.md src/*.c src/*.h src/Makefile tests/*.sh scripts/*.sh scripts/*.py examples/*.md examples/*.sam docs/*.md benchmarks/*.md CITATION.cff CHANGELOG.md Makefile .github/workflows/*.yml .gitignore'

if grep -n "$stale_name_upper\|$stale_name_lower\|$stale_version_prefix" $text_paths 2>/dev/null; then
  echo "Found stale pre-publication names or versions" >&2
  exit 1
fi

release_metadata_paths='README.md CITATION.cff docs/release_checklist.md'
for metadata_path in $release_metadata_paths; do
  if [ ! -s "$metadata_path" ]; then
    echo "Missing or empty release metadata file: $metadata_path" >&2
    exit 1
  fi
done
placeholder_pattern='github[.]com/your-org\|to be added after release\|update `CITATION[.]cff`\|add the Zenodo or institutional archive DOI'

if grep -n "$placeholder_pattern" $release_metadata_paths; then
  if [ "${REQUIRE_RELEASE_METADATA:-0}" = "1" ]; then
    echo "Found placeholder release metadata; replace unfinished repository or archive entries before release" >&2
    exit 1
  fi
  echo "WARNING: placeholder release metadata remains; run make check-release before tagging" >&2
fi

printf '%s\n' "Release metadata check passed for $expected_version"
