#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out_dir="$repo_dir/examples/output"
mkdir -p "$out_dir"

"$repo_dir/src/isoumi" \
  --bam "$repo_dir/examples/minimal.sam" \
  --out "$out_dir/minimal" \
  --tmp-dir "$out_dir/buckets" \
  --threads 1 \
  --buckets 4 \
  --ham 1 \
  --ratio 0.60 \
  --emit-tsv \
  --emit-explain

printf '%s\n' "Wrote example outputs to $out_dir"
