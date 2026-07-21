#!/bin/sh
set -eu

bin=${1:-src/isoumi}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-directional.XXXXXX")

cleanup() {
  rm -f "$tmpdir/input.sam"
  rm -f "$tmpdir/graph.sam"
  rm -f "$tmpdir/bad.stdout" "$tmpdir/bad.stderr"
  rm -f "$tmpdir/conflict.stdout" "$tmpdir/conflict.stderr"
  rm -f "$tmpdir/ratio.dedup.bam"
  rm -f "$tmpdir/ratio.molecules.tsv" "$tmpdir/ratio.assignments.tsv"
  rm -f "$tmpdir/ratio.corrections.tsv"
  rm -f "$tmpdir/directional.dedup.bam"
  rm -f "$tmpdir/directional.molecules.tsv" "$tmpdir/directional.assignments.tsv"
  rm -f "$tmpdir/directional.corrections.tsv"
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

cat > "$tmpdir/input.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
d001	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" --help 2>&1 | grep -q -- '--correction-method <ratio|directional>'

if "$bin" --bam "$tmpdir/input.sam" --out "$tmpdir/bad" \
    --correction-method unknown >"$tmpdir/bad.stdout" 2>"$tmpdir/bad.stderr"; then
  echo "unknown correction method should fail" >&2
  exit 1
fi
grep -q 'unknown --correction-method: unknown' "$tmpdir/bad.stderr"

if "$bin" --bam "$tmpdir/input.sam" --out "$tmpdir/conflict" \
    --min-merge-confidence 0.10 --correction-method directional \
    >"$tmpdir/conflict.stdout" 2>"$tmpdir/conflict.stderr"; then
  echo "directional with positive merge confidence should fail" >&2
  exit 1
fi
grep -q -- '--min-merge-confidence is incompatible with --correction-method directional' \
  "$tmpdir/conflict.stderr"

cat > "$tmpdir/graph.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
SAM

next_read=1
append_umi_reads() {
  cell=$1
  umi=$2
  count=$3
  i=1
  while [ "$i" -le "$count" ]; do
    printf 'g%03d\t0\tchr1\t101\t60\t10M\t*\t0\t0\tACGTACGTAA\tFFFFFFFFFF\tCB:Z:%s\tUR:Z:%s\tUY:Z:IIII\tGX:Z:GENE1\n' \
      "$next_read" "$cell" "$umi" >> "$tmpdir/graph.sam"
    next_read=$((next_read + 1))
    i=$((i + 1))
  done
}

append_umi_reads CELL_CHAIN AAAA 5
append_umi_reads CELL_CHAIN AAAT 3
append_umi_reads CELL_CHAIN AATT 2
append_umi_reads CELL_CHAIN ATTT 2
append_umi_reads CELL_BOUNDARY CCCC 3
append_umi_reads CELL_BOUNDARY CCCT 2
append_umi_reads CELL_TIE GGGA 1
append_umi_reads CELL_TIE GGGG 1
append_umi_reads CELL_DIRECT TTTT 10
append_umi_reads CELL_DIRECT TTTA 1

"$bin" \
  --bam "$tmpdir/graph.sam" \
  --out "$tmpdir/ratio" \
  --threads 1 \
  --buckets 2 \
  --ham 1 \
  --emit-tsv \
  --emit-explain >/dev/null

"$bin" \
  --bam "$tmpdir/graph.sam" \
  --out "$tmpdir/directional" \
  --threads 1 \
  --buckets 2 \
  --ham 1 \
  --correction-method directional \
  --emit-tsv \
  --emit-explain >/dev/null

awk -F '\t' 'NR>1 && $2=="AATT" && $3=="AATT" \
  {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/ratio.corrections.tsv"

expected_header=$(printf 'key\traw_umi\tcorr_umi\traw_count\tseed_count\thamming\traw_avgq\tseed_avgq\tconfidence\treason\tbucket\tmethod\tparent_umi\tparent_count\tedge_hamming\tpath_length')
test "$(sed -n '1p' "$tmpdir/directional.corrections.tsv")" = "$expected_header"

awk -F '\t' 'NR>1 && $2=="AATT" && $3=="AAAA" && $6==2 && \
  $10=="directional" && $12=="directional" && $13=="AAAT" && \
  $14==3 && $15==1 && $16==2 {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/directional.corrections.tsv"

awk -F '\t' 'NR>1 && $2=="ATTT" && $3=="ATTT" && $16==0 \
  {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/directional.corrections.tsv"

awk -F '\t' 'NR>1 && $2=="CCCT" && $3=="CCCC" && $16==1 \
  {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/directional.corrections.tsv"

awk -F '\t' 'NR>1 && $2=="GGGG" && $3=="GGGA" \
  {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/directional.corrections.tsv"

awk -F '\t' 'NR>1 && $2=="TTTA" && $3=="TTTT" && $12=="ratio" && \
  $13=="TTTT" && $14==10 && $15==1 && $16==1 \
  {found=1} END {exit found ? 0 : 1}' \
  "$tmpdir/ratio.corrections.tsv"
