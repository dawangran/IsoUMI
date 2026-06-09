#!/bin/sh
set -eu

bin=${1:-src/isoumi}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-smoke.XXXXXX")

cleanup() {
  rm -f "$tmpdir/smoke.sam" "$tmpdir/tie.sam" \
        "$tmpdir/out.dedup.bam" \
        "$tmpdir/out.molecules.tsv" "$tmpdir/out.assignments.tsv" \
        "$tmpdir/out.corrections.tsv" \
        "$tmpdir/tie.dedup.bam" \
        "$tmpdir/tie.molecules.tsv" "$tmpdir/tie.assignments.tsv" \
        "$tmpdir/tie.corrections.tsv"
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

cat > "$tmpdir/smoke.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r001	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
r002	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
r003	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:!!!!	GX:Z:GENE1
r004	0	chr1	301	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:TTTT	UY:Z:IIII	GX:Z:GENE1
r005	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	UR:Z:AAAG	UY:Z:!!!!	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/smoke.sam" \
  --out "$tmpdir/out" \
  --tmp-dir "$tmpdir/buckets" \
  --threads 1 \
  --buckets 4 \
  --ham 1 \
  --ratio 0.60 \
  --mol-tag MI \
  --emit-tsv \
  --emit-explain >/dev/null

test -s "$tmpdir/out.dedup.bam"
test -s "$tmpdir/out.molecules.tsv"
test -s "$tmpdir/out.assignments.tsv"
test -s "$tmpdir/out.corrections.tsv"
grep -q 'AAAT	AAAA' "$tmpdir/out.corrections.tsv"
grep -q 'count+quality' "$tmpdir/out.corrections.tsv"
grep -q '^r001	CELL1|' "$tmpdir/out.assignments.tsv"
grep -q '^r002	CELL1|' "$tmpdir/out.assignments.tsv"
grep -q '^r003	CELL1|' "$tmpdir/out.assignments.tsv"
grep -q '^r004	CELL1|' "$tmpdir/out.assignments.tsv"
if grep -Eq 'CB=NA|^r005	' "$tmpdir/out.corrections.tsv" "$tmpdir/out.assignments.tsv"; then
  echo "missing-cell read should pass through without correction reports" >&2
  exit 1
fi

cat > "$tmpdir/tie.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
t001	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:!!!!	GX:Z:GENE1
t002	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/tie.sam" \
  --out "$tmpdir/tie" \
  --tmp-dir "$tmpdir/tie_buckets" \
  --threads 1 \
  --buckets 2 \
  --ham 1 \
  --ratio 1.00 \
  --no-quality-aware \
  --emit-tsv \
  --emit-explain >/dev/null

grep -q 'AAAT	AAAA' "$tmpdir/tie.corrections.tsv"
