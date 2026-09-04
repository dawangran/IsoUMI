#!/bin/sh
set -eu

bin=${1:-src/isoumi}

if ! command -v samtools >/dev/null 2>&1; then
  echo "SKIP: samtools not found; non-primary alignment test not run" >&2
  exit 0
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-non-primary.XXXXXX")

cleanup() {
  rm -f "$tmpdir/input.sam"
  rm -f "$tmpdir/out.dedup.bam"
  rm -f "$tmpdir/out.molecules.tsv" "$tmpdir/out.assignments.tsv"
  rm -f "$tmpdir/out.corrections.tsv" "$tmpdir/out.view"
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

cat > "$tmpdir/input.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
p1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
p1	2048	chr1	201	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
p1	2048	chr1	301	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
p1	2048	chr1	401	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
p1	2048	chr1	501	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
p2	0	chr1	101	50	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
p3	0	chr1	101	40	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
orphan	2048	chr1	701	30	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:CCCC	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/input.sam" \
  --out "$tmpdir/out" \
  --threads 2 \
  --buckets 4 \
  --ham 1 \
  --ratio 0.60 \
  --mol-tag MI \
  --emit-tsv \
  --emit-explain >/dev/null

samtools view "$tmpdir/out.dedup.bam" > "$tmpdir/out.view"

awk -F '\t' '$1 == "p1" {
  n++; ub=0; da=0; mi=0
  for (i=12; i<=NF; i++) {
    if ($i == "UB:Z:AAAT") ub=1
    if ($i ~ /^DA:[A-Za-z]:0$/) da=1
    if ($i ~ /^MI:Z:/) mi=1
  }
  if (!ub || !da || !mi) bad=1
} END {exit n == 5 && !bad ? 0 : 1}' "$tmpdir/out.view"

awk -F '\t' '$1 == "p2" || $1 == "p3" {
  n++; ub=0; da=0
  for (i=12; i<=NF; i++) {
    if ($i == "UB:Z:AAAT") ub=1
    if ($i ~ /^DA:[A-Za-z]:1$/) da=1
  }
  if (!ub || !da || $2 != 1024) bad=1
} END {exit n == 2 && !bad ? 0 : 1}' "$tmpdir/out.view"

awk -F '\t' '$1 == "orphan" {
  found=1; ub=0; da=0
  for (i=12; i<=NF; i++) {
    if ($i == "UB:Z:CCCC") ub=1
    if ($i ~ /^DA:[A-Za-z]:0$/) da=1
  }
  if (!ub || !da) exit 2
} END {exit found ? 0 : 1}' "$tmpdir/out.view"

test "$(awk 'NR > 1 {n++} END {print n + 0}' "$tmpdir/out.assignments.tsv")" = "3"
test "$(awk 'NR > 1 {n++} END {print n + 0}' "$tmpdir/out.molecules.tsv")" = "1"
awk -F '\t' 'NR > 1 && $2 == "AAAT" && $3 == 3 {found=1}
  END {exit found ? 0 : 1}' "$tmpdir/out.molecules.tsv"
awk -F '\t' 'NR > 1 && $2 == "AAAA" && $3 == "AAAT" && $4 == 1 && $5 == 2 {found=1}
  END {exit found ? 0 : 1}' "$tmpdir/out.corrections.tsv"

printf '%s\n' "Non-primary alignment test passed"
