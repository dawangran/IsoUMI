#!/bin/sh
set -eu

version=$(cat VERSION)
expected_name="IsoUMI"
expected_version="$expected_name $version"

grep -q "$expected_version" src/cli.c
grep -q "$expected_version" README.md
grep -q "version: \"$version\"" CITATION.cff
grep -q "## $expected_name $version" CHANGELOG.md
grep -q "BIN = isoumi" src/Makefile
grep -q "src/isoumi" .gitignore

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

printf '%s\n' "Release metadata check passed for $expected_version"
