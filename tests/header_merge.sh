#!/bin/sh
set -eu

bin=${1:-src/isoumi}

if ! command -v samtools >/dev/null 2>&1; then
  echo "SKIP: samtools not found; header merge test not run" >&2
  exit 0
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-header-merge.XXXXXX")

cleanup() {
  rm -f "$tmpdir/input_a.sam"
  rm -f "$tmpdir/input_b.sam"
  rm -f "$tmpdir/rg_a.sam"
  rm -f "$tmpdir/rg_b.sam"
  rm -f "$tmpdir/ref_a.sam"
  rm -f "$tmpdir/ref_b.sam"
  rm -f "$tmpdir/out.dedup.bam"
  rm -f "$tmpdir/out.stderr"
  rm -f "$tmpdir/out.header"
  rm -f "$tmpdir/out.view"
  rm -f "$tmpdir/strip.dedup.bam"
  rm -f "$tmpdir/strip.header"
  rm -f "$tmpdir/strip.view"
  rm -f "$tmpdir/strip.bucket.header"
  rm -f "$tmpdir/strip.bucket.view"
  rm -f "$tmpdir/rg.stderr"
  rm -f "$tmpdir/ref.stderr"
  rm -f "$tmpdir/buckets/bucket_00000.bam"
  rm -f "$tmpdir/buckets/bucket_00000.dedup.bam"
  rm -f "$tmpdir/strip_buckets/bucket_00000.bam"
  rm -f "$tmpdir/strip_buckets/bucket_00000.dedup.bam"
  rmdir "$tmpdir/buckets" 2>/dev/null || true
  rmdir "$tmpdir/strip_buckets" 2>/dev/null || true
  rmdir "$tmpdir/rg_buckets" 2>/dev/null || true
  rmdir "$tmpdir/ref_buckets" 2>/dev/null || true
  rmdir "$tmpdir" 2>/dev/null || true
}
trap cleanup EXIT HUP INT TERM

cat > "$tmpdir/input_a.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@PG	ID:minimap2	PN:minimap2	VN:2.26	CL:minimap2 -ax splice ref.fa reads_a.fq
@PG	ID:samtools	PN:samtools	PP:minimap2	VN:1.20	CL:samtools sort input_a.sam
a1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	PG:Z:samtools
SAM

cat > "$tmpdir/input_b.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@PG	ID:minimap2	PN:minimap2	VN:2.28	CL:minimap2 -ax splice ref.fa reads_b.fq
@PG	ID:samtools	PN:samtools	PP:minimap2	VN:1.21	CL:samtools sort input_b.sam
b1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	PG:Z:samtools
SAM

"$bin" \
  --bam "$tmpdir/input_a.sam" \
  --bam "$tmpdir/input_b.sam" \
  --out "$tmpdir/out" \
  --tmp-dir "$tmpdir/buckets" \
  --threads 1 \
  --buckets 1 \
  --no-gene \
  --no-quality-aware \
  2> "$tmpdir/out.stderr"

test "$(samtools view -c "$tmpdir/out.dedup.bam")" = "2"
samtools view --no-PG -H "$tmpdir/out.dedup.bam" > "$tmpdir/out.header"
samtools view "$tmpdir/out.dedup.bam" > "$tmpdir/out.view"
test "$(awk '$1 == "@PG" {n++} END {print n + 0}' "$tmpdir/out.header")" = "5"
awk -F '\t' '$1 == "@PG" {
  id = ""; vn = ""
  for (i = 2; i <= NF; i++) {
    if ($i ~ /^ID:/) id = substr($i, 4)
    if ($i ~ /^VN:/) vn = substr($i, 4)
  }
  if (id == "minimap2" && vn == "2.26") first = 1
  if (id == "minimap2.isoumi.input000001" && vn == "2.28") second = 1
} END {exit first && second ? 0 : 1}' "$tmpdir/out.header"
awk -F '\t' '$1 == "@PG" {
  id = ""; pn = ""; vn = ""; cl = ""
  for (i = 2; i <= NF; i++) {
    if ($i ~ /^ID:/) id = substr($i, 4)
    if ($i ~ /^PN:/) pn = substr($i, 4)
    if ($i ~ /^VN:/) vn = substr($i, 4)
    if ($i ~ /^CL:/) cl = substr($i, 4)
  }
  if (id == "IsoUMI" && pn == "IsoUMI" && vn == "0.1.1" && cl != "") found = 1
} END {exit found ? 0 : 1}' "$tmpdir/out.header"
awk -F '\t' '$1 == "@PG" {
  id = ""; pp = ""
  for (i = 2; i <= NF; i++) {
    if ($i ~ /^ID:/) id = substr($i, 4)
    if ($i ~ /^PP:/) pp = substr($i, 4)
  }
  if (id == "samtools.isoumi.input000001" && pp == "minimap2.isoumi.input000001") found = 1
} END {exit found ? 0 : 1}' "$tmpdir/out.header"
awk -F '\t' '$1 == "a1" {for (i=12; i<=NF; i++) if ($i == "PG:Z:samtools") found=1}
  END {exit found ? 0 : 1}' "$tmpdir/out.view"
awk -F '\t' '$1 == "b1" {for (i=12; i<=NF; i++) if ($i == "PG:Z:samtools.isoumi.input000001") found=1}
  END {exit found ? 0 : 1}' "$tmpdir/out.view"
grep -q 'renaming @PG ID=minimap2 from input 1 to minimap2.isoumi.input000001' "$tmpdir/out.stderr"

"$bin" \
  --bam "$tmpdir/input_a.sam" \
  --bam "$tmpdir/input_b.sam" \
  --out "$tmpdir/strip" \
  --tmp-dir "$tmpdir/strip_buckets" \
  --threads 1 \
  --buckets 1 \
  --no-gene \
  --no-quality-aware \
  --strip-pg \
  --keep-tmp \
  2>/dev/null

test "$(samtools view -c "$tmpdir/strip.dedup.bam")" = "2"
samtools view --no-PG -H "$tmpdir/strip.dedup.bam" > "$tmpdir/strip.header"
samtools view "$tmpdir/strip.dedup.bam" > "$tmpdir/strip.view"
samtools view --no-PG -H "$tmpdir/strip_buckets/bucket_00000.bam" > "$tmpdir/strip.bucket.header"
samtools view "$tmpdir/strip_buckets/bucket_00000.bam" > "$tmpdir/strip.bucket.view"
if awk '$1 == "@PG" {found=1} END {exit found ? 0 : 1}' "$tmpdir/strip.header"; then
  echo "--strip-pg left an @PG header line" >&2
  exit 1
fi
if awk -F '\t' '{for (i=12; i<=NF; i++) if ($i ~ /^PG:/) found=1}
  END {exit found ? 0 : 1}' "$tmpdir/strip.view"; then
  echo "--strip-pg left a record-level PG tag" >&2
  exit 1
fi
if awk '$1 == "@PG" {found=1} END {exit found ? 0 : 1}' "$tmpdir/strip.bucket.header" ||
   awk -F '\t' '{for (i=12; i<=NF; i++) if ($i ~ /^PG:/) found=1}
     END {exit found ? 0 : 1}' "$tmpdir/strip.bucket.view"; then
  echo "--strip-pg left PG provenance in an intermediate bucket BAM" >&2
  exit 1
fi

cat > "$tmpdir/rg_a.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@RG	ID:group1	SM:sample_a
ra1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	RG:Z:group1
SAM

cat > "$tmpdir/rg_b.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@RG	ID:group1	SM:sample_b
rb1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	RG:Z:group1
SAM


if "$bin" \
  --bam "$tmpdir/rg_a.sam" \
  --bam "$tmpdir/rg_b.sam" \
  --out "$tmpdir/rg" \
  --tmp-dir "$tmpdir/rg_buckets" \
  --threads 1 \
  --buckets 1 \
  2> "$tmpdir/rg.stderr"; then
  echo "conflicting @RG definitions should remain fatal" >&2
  exit 1
fi
grep -q 'conflicting @RG header line for ID=group1' "$tmpdir/rg.stderr"

cat > "$tmpdir/ref_a.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000	M5:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
ma1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA
SAM

cat > "$tmpdir/ref_b.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000	M5:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
mb1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA
SAM

if "$bin" \
  --bam "$tmpdir/ref_a.sam" \
  --bam "$tmpdir/ref_b.sam" \
  --out "$tmpdir/ref" \
  --tmp-dir "$tmpdir/ref_buckets" \
  --threads 1 \
  --buckets 1 \
  2> "$tmpdir/ref.stderr"; then
  echo "different @SQ M5 values should be rejected" >&2
  exit 1
fi
grep -q 'reference metadata mismatch for chr1: @SQ M5:' "$tmpdir/ref.stderr"

printf '%s\n' "Header merge test passed"
