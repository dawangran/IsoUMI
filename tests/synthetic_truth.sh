#!/bin/sh
set -eu

bin=${1:-src/isoumi}
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-synthetic.XXXXXX")

remove_file() {
  rm -f "$1"
}

remove_case_outputs() {
  prefix="$tmpdir/$1"
  remove_file "$prefix.dedup.bam"
  remove_file "$prefix.molecules.tsv"
  remove_file "$prefix.assignments.tsv"
  remove_file "$prefix.corrections.tsv"
  remove_file "$prefix.stdout.log"
  remove_file "$prefix.stderr.log"
  remove_file "$prefix.time.txt"
}

cleanup() {
  remove_file "$tmpdir/synthetic_truth.sam"
  remove_file "$tmpdir/synthetic_truth.truth.tsv"
  remove_case_outputs default
  remove_case_outputs baseline
  remove_case_outputs endbin50
  remove_case_outputs conservative
  remove_case_outputs baseline_no_structure
  remove_file "$tmpdir/metrics.tsv"
  remove_file "$tmpdir/summary.tsv"
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

python3 "$repo_dir/scripts/generate_synthetic_truth.py" \
  --out-dir "$tmpdir" \
  --prefix synthetic_truth >/dev/null

test -s "$tmpdir/synthetic_truth.sam"
test -s "$tmpdir/synthetic_truth.truth.tsv"
grep -q '^qname	scenario	true_molecule' "$tmpdir/synthetic_truth.truth.tsv"
grep -q 'umi_error_same_junction' "$tmpdir/synthetic_truth.truth.tsv"
grep -q 'same_umi_different_junction' "$tmpdir/synthetic_truth.truth.tsv"
grep -q 'missing_cell_passthrough' "$tmpdir/synthetic_truth.truth.tsv"

"$bin" \
  --bam "$tmpdir/synthetic_truth.sam" \
  --out "$tmpdir/default" \
  --threads 1 \
  --buckets 8 \
  --ham 1 \
  --ratio 0.60 \
  --emit-tsv \
  --emit-explain >/dev/null

"$bin" \
  --bam "$tmpdir/synthetic_truth.sam" \
  --out "$tmpdir/baseline" \
  --threads 1 \
  --buckets 8 \
  --ham 1 \
  --ratio 0.60 \
  --no-structure \
  --emit-tsv \
  --emit-explain >/dev/null

test -s "$tmpdir/default.dedup.bam"
test -s "$tmpdir/default.molecules.tsv"
test -s "$tmpdir/default.assignments.tsv"
test -s "$tmpdir/default.corrections.tsv"
test -s "$tmpdir/baseline.molecules.tsv"

grep -q 'AAAT	AAAA' "$tmpdir/default.corrections.tsv"
grep -q 'TTTA	TTTT' "$tmpdir/default.corrections.tsv"
if grep -Eq 'missing_cell|MOL_NO_CELL|CB=NA' "$tmpdir/default.corrections.tsv" "$tmpdir/default.assignments.tsv"; then
  echo "synthetic missing-cell reads should not appear in correction or assignment reports" >&2
  exit 1
fi

default_molecules=$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/default.molecules.tsv")
baseline_molecules=$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/baseline.molecules.tsv")

test "$default_molecules" = "6"
test "$baseline_molecules" = "5"

if ! awk -F '\t' 'NR>1 && $1 ~ /SJ=/ && $2 == "CCCC" {n++} END {exit n == 2 ? 0 : 1}' "$tmpdir/default.molecules.tsv"; then
  echo "default synthetic run should keep same-UMI different-junction molecules separate" >&2
  exit 1
fi

if ! awk -F '\t' 'NR>1 && $1 ~ /CTX=NA/ && $2 == "CCCC" && $3 == 6 {found=1} END {exit found ? 0 : 1}' "$tmpdir/baseline.molecules.tsv"; then
  echo "no-structure baseline should collapse same-UMI different-junction molecules" >&2
  exit 1
fi

sh "$repo_dir/scripts/benchmark_isoumi.sh" \
  --bin "$bin" \
  --dataset synthetic_truth \
  --threads 1 \
  --buckets 8 \
  "$tmpdir/synthetic_truth.sam" \
  "$tmpdir" >/dev/null

test -s "$tmpdir/metrics.tsv"
grep -q '^dataset	label	version	input_path	input_format	input_size_bytes	input_reads	input_cell_barcodes	input_umis' "$tmpdir/metrics.tsv"
grep -q '^synthetic_truth	default	' "$tmpdir/metrics.tsv"
grep -q '^synthetic_truth	baseline_no_structure	' "$tmpdir/metrics.tsv"

python3 "$repo_dir/scripts/summarize_benchmark_reports.py" "$tmpdir" >/dev/null
test -s "$tmpdir/summary.tsv"
grep -q '^synthetic_truth	default	' "$tmpdir/summary.tsv"
grep -q '^synthetic_truth	baseline_no_structure	' "$tmpdir/summary.tsv"
awk -F '\t' 'NR==1 {
  for (i=1; i<=NF; i++) h[$i]=i
  exit !(h["dataset"] && h["input_reads"] && h["assigned_cell_barcodes"] && h["peak_rss_mb"] && h["correction_rate"])
}' "$tmpdir/summary.tsv"
