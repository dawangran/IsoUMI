#!/bin/sh
set -eu

if [ $# -lt 2 ]; then
  echo "Usage: sh scripts/benchmark_isoumi.sh <input.bam|input.sam> <out-dir>" >&2
  exit 2
fi

input=$1
out_dir=$2
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
bin="$repo_dir/src/isoumi"
metrics="$out_dir/metrics.tsv"

if [ ! -x "$bin" ]; then
  echo "IsoUMI binary not found; run make first" >&2
  exit 2
fi

mkdir -p "$out_dir"
printf 'label\tversion\telapsed_seconds\tcommand\n' > "$metrics"

run_case() {
  label=$1
  shift
  prefix="$out_dir/$label"
  tmp="$out_dir/${label}.tmp"
  start=$(date +%s)
  "$bin" --bam "$input" --out "$prefix" --tmp-dir "$tmp" --threads 4 --buckets 64 "$@"
  end=$(date +%s)
  elapsed=$((end - start))
  version=$("$bin" --version 2>&1)
  printf '%s\t%s\t%s\t%s\n' "$label" "$version" "$elapsed" "$bin --bam $input --out $prefix --threads 4 --buckets 64 $*" >> "$metrics"
}

run_case default --emit-tsv --emit-explain
run_case endbin50 --end-bin 50 --emit-tsv --emit-explain
run_case conservative --ham 1 --ratio 0.05 --min-merge-confidence 0.60 --emit-tsv --emit-explain

printf '%s\n' "Wrote benchmark metrics to $metrics"
