#!/bin/sh
set -eu

bin=${1:-src/isoumi}

if ! command -v samtools >/dev/null 2>&1; then
  echo "SKIP: samtools not found; BAM tag semantics test not run" >&2
  exit 0
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-bam-tags.XXXXXX")

remove_file() {
  rm -f "$1"
}

cleanup() {
  remove_file "$tmpdir/tags.sam"
  remove_file "$tmpdir/tags.view"
  remove_file "$tmpdir/tags.dedup.bam"
  remove_file "$tmpdir/no_flag.sam"
  remove_file "$tmpdir/no_flag.view"
  remove_file "$tmpdir/no_flag.dedup.bam"
  remove_file "$tmpdir/input_a.sam"
  remove_file "$tmpdir/input_b.sam"
  remove_file "$tmpdir/isolated.view"
  remove_file "$tmpdir/isolated.dedup.bam"
  remove_file "$tmpdir/isolated.molecules.tsv"
  remove_file "$tmpdir/isolated.assignments.tsv"
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

assert_read() {
  view=$1
  qname=$2
  expected_flag=$3
  expected_da=$4
  expected_ub=$5
  awk -v q="$qname" -v flag="$expected_flag" -v da="$expected_da" -v ub="$expected_ub" '
    $1 == q {
      found = 1
      if ($2 != flag) exit 2
      has_da = 0
      has_ub = 0
      has_zi = 0
      for (i = 12; i <= NF; i++) {
        if ($i ~ "^DA:[A-Za-z]:" da "$") has_da = 1
        if ($i == "UB:Z:" ub) has_ub = 1
        if ($i ~ "^zi:Z:") has_zi = 1
      }
      if (!has_da) exit 3
      if (!has_ub) exit 4
      if (has_zi) exit 5
    }
    END { if (!found) exit 1 }
  ' "$view"
}

cat > "$tmpdir/tags.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r001	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
r002	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
r003	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:!!!!	GX:Z:GENE1
r004	0	chr1	301	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:TTTT	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/tags.sam" \
  --out "$tmpdir/tags" \
  --threads 1 \
  --buckets 4 \
  --ham 1 \
  --ratio 0.60 >/dev/null

samtools view "$tmpdir/tags.dedup.bam" > "$tmpdir/tags.view"
assert_read "$tmpdir/tags.view" r001 0 0 AAAA
assert_read "$tmpdir/tags.view" r002 1024 1 AAAA
assert_read "$tmpdir/tags.view" r003 1024 1 AAAA
assert_read "$tmpdir/tags.view" r004 0 0 TTTT

cat > "$tmpdir/no_flag.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
keep1	1024	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:CCCC	UY:Z:IIII	GX:Z:GENE1
keep2	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:CCCC	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/no_flag.sam" \
  --out "$tmpdir/no_flag" \
  --threads 1 \
  --buckets 2 \
  --no-bam-dup-flag >/dev/null

samtools view "$tmpdir/no_flag.dedup.bam" > "$tmpdir/no_flag.view"
assert_read "$tmpdir/no_flag.view" keep1 1024 0 CCCC
assert_read "$tmpdir/no_flag.view" keep2 0 1 CCCC

cat > "$tmpdir/input_a.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
a1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELLX	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE1
SAM

cat > "$tmpdir/input_b.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
b1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELLX	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/input_a.sam" \
  --bam "$tmpdir/input_b.sam" \
  --out "$tmpdir/isolated" \
  --threads 1 \
  --buckets 2 \
  --isolate-inputs \
  --emit-tsv >/dev/null

samtools view "$tmpdir/isolated.dedup.bam" > "$tmpdir/isolated.view"
if awk '{ for (i = 12; i <= NF; i++) if ($i ~ /^zi:Z:/) found = 1 } END { exit found ? 0 : 1 }' "$tmpdir/isolated.view"; then
  echo "internal input-scope tag zi leaked into final BAM" >&2
  exit 1
fi
grep -q 'IN=input000000' "$tmpdir/isolated.molecules.tsv"
grep -q 'IN=input000001' "$tmpdir/isolated.molecules.tsv"
