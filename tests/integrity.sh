#!/bin/sh
set -eu

bin=${1:-src/isoumi}

if ! command -v samtools >/dev/null 2>&1; then
  echo "SKIP: samtools not found; integrity regression test not run" >&2
  exit 0
fi

tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-integrity.XXXXXX")
cleanup() {
  rm -rf "$tmpdir"
}
trap cleanup EXIT HUP INT TERM

cat > "$tmpdir/base.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
base	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	GX:Z:GENE1
SAM

samtools view -bS -o "$tmpdir/alias.dedup.bam" "$tmpdir/base.sam"
alias_before=$(cksum < "$tmpdir/alias.dedup.bam")
if "$bin" --bam "$tmpdir/alias.dedup.bam" --out "$tmpdir/alias" \
    --tmp-dir "$tmpdir/alias_buckets" --threads 1 --buckets 1 \
    >"$tmpdir/alias.stdout" 2>"$tmpdir/alias.stderr"; then
  echo "input/output aliases should be rejected" >&2
  exit 1
fi
test "$(cksum < "$tmpdir/alias.dedup.bam")" = "$alias_before"
test ! -e "$tmpdir/alias_buckets"
grep -q 'refusing to overwrite input' "$tmpdir/alias.stderr"

cat > "$tmpdir/keys.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
k1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:A|GX=B	UR:Z:AAAA	GX:Z:C
k2	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:A	UR:Z:AAAA	GX:Z:B|GX=C
stale	1024	chr1	301	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	GX:Z:GENE1	UB:Z:STALE	MI:Z:OLD
SAM

"$bin" --bam "$tmpdir/keys.sam" --out "$tmpdir/keys" \
  --tmp-dir "$tmpdir/key_buckets" --threads 1 --buckets 2 \
  --mol-tag MI --emit-tsv >/dev/null
samtools view "$tmpdir/keys.dedup.bam" > "$tmpdir/keys.view"
test "$(awk 'NR > 1 {n++} END {print n + 0}' "$tmpdir/keys.molecules.tsv")" = "2"
grep -q 'CB=A%7CGX%3DB|GX=C' "$tmpdir/keys.molecules.tsv"
grep -q 'CB=A|GX=B%7CGX%3DC' "$tmpdir/keys.molecules.tsv"
awk -F '\t' '$1 == "k1" || $1 == "k2" {
  n++
  if ($2 != 0) bad=1
  da=0
  for (i=12; i<=NF; i++) if ($i ~ /^DA:[A-Za-z]:0$/) da=1
  if (!da) bad=1
} END {exit n == 2 && !bad ? 0 : 1}' "$tmpdir/keys.view"
awk -F '\t' '$1 == "stale" {
  found=1
  if ($2 != 0) bad=1
  da=0
  for (i=12; i<=NF; i++) {
    if ($i ~ /^(UB|MI):/) bad=1
    if ($i ~ /^DA:[A-Za-z]:0$/) da=1
  }
  if (!da) bad=1
} END {exit found && !bad ? 0 : 1}' "$tmpdir/keys.view"

if "$bin" --bam "$tmpdir/base.sam" --out "$tmpdir/reserved" \
    --strip-pg --umi-out PG >"$tmpdir/reserved.stdout" 2>"$tmpdir/reserved.stderr"; then
  echo "PG should be reserved for SAM provenance" >&2
  exit 1
fi
grep -q 'reserved SAM provenance tag PG' "$tmpdir/reserved.stderr"

cat > "$tmpdir/ref_none.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
SAM
cat > "$tmpdir/ref_a.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000	M5:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa	AS:buildA
SAM
cat > "$tmpdir/ref_b.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000	M5:bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb	AS:buildB
SAM

"$bin" --bam "$tmpdir/ref_none.sam" --bam "$tmpdir/ref_a.sam" \
  --out "$tmpdir/ref_ok" --tmp-dir "$tmpdir/ref_ok_buckets" \
  --threads 1 --buckets 1 >/dev/null
samtools view --no-PG -H "$tmpdir/ref_ok.dedup.bam" > "$tmpdir/ref_ok.header"
awk -F '\t' '$1 == "@SQ" {
  for (i=2; i<=NF; i++) {
    if ($i == "M5:aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa") m5=1
    if ($i == "AS:buildA") as=1
  }
} END {exit m5 && as ? 0 : 1}' "$tmpdir/ref_ok.header"

if "$bin" --bam "$tmpdir/ref_none.sam" --bam "$tmpdir/ref_a.sam" \
    --bam "$tmpdir/ref_b.sam" --out "$tmpdir/ref_bad" \
    --tmp-dir "$tmpdir/ref_bad_buckets" --threads 1 --buckets 1 \
    >"$tmpdir/ref_bad.stdout" 2>"$tmpdir/ref_bad.stderr"; then
  echo "three-way @SQ conflicts should be rejected" >&2
  exit 1
fi
test ! -e "$tmpdir/ref_bad_buckets"
grep -q 'reference metadata mismatch for chr1: @SQ M5:' "$tmpdir/ref_bad.stderr"

cat > "$tmpdir/pg_root.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
@PG	ID:root	PN:root
SAM
cat > "$tmpdir/pg_child.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
@PG	ID:child	PN:child	PP:root
SAM
if "$bin" --bam "$tmpdir/pg_root.sam" --bam "$tmpdir/pg_child.sam" \
    --out "$tmpdir/pg_bad" --tmp-dir "$tmpdir/pg_bad_buckets" \
    --threads 1 --buckets 1 >"$tmpdir/pg_bad.stdout" 2>"$tmpdir/pg_bad.stderr"; then
  echo "cross-input @PG parent links should be rejected" >&2
  exit 1
fi
test ! -e "$tmpdir/pg_bad_buckets"
grep -q 'references missing parent PP=root' "$tmpdir/pg_bad.stderr"

cat > "$tmpdir/unmapped_supp.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
same	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	GX:Z:GENE1
same	3076	*	0	0	*	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	GX:Z:GENE1	UB:Z:STALE	MI:Z:OLD
SAM
"$bin" --bam "$tmpdir/unmapped_supp.sam" --out "$tmpdir/unmapped_supp" \
  --tmp-dir "$tmpdir/unmapped_supp_buckets" --threads 1 --buckets 1 \
  --mol-tag MI >/dev/null
samtools view "$tmpdir/unmapped_supp.dedup.bam" > "$tmpdir/unmapped_supp.view"
awk -F '\t' '$1 == "same" && $2 == 2052 {
  found=1
  ub=0
  da=0
  for (i=12; i<=NF; i++) {
    if ($i == "UB:Z:AAAA") ub=1
    if ($i ~ /^DA:[A-Za-z]:0$/) da=1
    if ($i ~ /^MI:/) bad=1
  }
  if (!ub || !da) bad=1
} END {exit found && !bad ? 0 : 1}' "$tmpdir/unmapped_supp.view"

cat > "$tmpdir/missing_cb_nonprimary.sam" <<'SAM'
@HD	VN:1.6	SO:unknown
@SQ	SN:chr1	LN:2000
q1	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	GX:Z:GENE1
q1	2048	chr1	201	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	UR:Z:AAAA	GX:Z:GENE1
SAM
if "$bin" --bam "$tmpdir/missing_cb_nonprimary.sam" --out "$tmpdir/missing_cb" \
    --tmp-dir "$tmpdir/missing_cb_buckets" --threads 1 --buckets 2 \
    >"$tmpdir/missing_cb.stdout" 2>"$tmpdir/missing_cb.stderr"; then
  echo "mapped non-primary records without CB should fail" >&2
  exit 1
fi
test ! -e "$tmpdir/missing_cb_buckets"
grep -q 'lacks CB; cannot reliably inherit' "$tmpdir/missing_cb.stderr"

if find "$tmpdir" -name '*.partial.*' -print | grep -q .; then
  echo "staging files leaked after completed or failed runs" >&2
  exit 1
fi

printf '%s\n' "Integrity regression test passed"
