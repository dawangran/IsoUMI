#!/bin/sh
set -eu

bin=${1:-src/isoumi}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-smoke.XXXXXX")

remove_file() {
  rm -f "$1"
}

cleanup() {
  remove_file "$tmpdir/smoke.sam"
  remove_file "$tmpdir/tie.sam"
  remove_file "$tmpdir/chain.sam"
  remove_file "$tmpdir/jitter.sam"
  remove_file "$tmpdir/structure.sam"
  remove_file "$tmpdir/rep.sam"
  remove_file "$tmpdir/endbin.sam"
  remove_file "$tmpdir/multi_a.sam"
  remove_file "$tmpdir/multi_b.sam"
  remove_file "$tmpdir/out.dedup.bam"
  remove_file "$tmpdir/out.molecules.tsv"
  remove_file "$tmpdir/out.assignments.tsv"
  remove_file "$tmpdir/out.corrections.tsv"
  remove_file "$tmpdir/tie.dedup.bam"
  remove_file "$tmpdir/tie.molecules.tsv"
  remove_file "$tmpdir/tie.assignments.tsv"
  remove_file "$tmpdir/tie.corrections.tsv"
  remove_file "$tmpdir/chain.dedup.bam"
  remove_file "$tmpdir/chain.corrections.tsv"
  remove_file "$tmpdir/jitter.dedup.bam"
  remove_file "$tmpdir/jitter.corrections.tsv"
  remove_file "$tmpdir/structure_default.dedup.bam"
  remove_file "$tmpdir/structure_default.molecules.tsv"
  remove_file "$tmpdir/structure_default.assignments.tsv"
  remove_file "$tmpdir/structure_baseline.dedup.bam"
  remove_file "$tmpdir/structure_baseline.molecules.tsv"
  remove_file "$tmpdir/structure_baseline.assignments.tsv"
  remove_file "$tmpdir/rep.dedup.bam"
  remove_file "$tmpdir/rep.molecules.tsv"
  remove_file "$tmpdir/rep.assignments.tsv"
  remove_file "$tmpdir/endbin_off.dedup.bam"
  remove_file "$tmpdir/endbin_off.molecules.tsv"
  remove_file "$tmpdir/endbin_off.assignments.tsv"
  remove_file "$tmpdir/endbin_on.dedup.bam"
  remove_file "$tmpdir/endbin_on.molecules.tsv"
  remove_file "$tmpdir/endbin_on.assignments.tsv"
  remove_file "$tmpdir/multi_shared.dedup.bam"
  remove_file "$tmpdir/multi_shared.molecules.tsv"
  remove_file "$tmpdir/multi_shared.assignments.tsv"
  remove_file "$tmpdir/multi_isolated.dedup.bam"
  remove_file "$tmpdir/multi_isolated.molecules.tsv"
  remove_file "$tmpdir/multi_isolated.assignments.tsv"
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

cat > "$tmpdir/chain.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
c001	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
c002	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
c003	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
c004	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
c005	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
c006	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
c007	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
c008	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:AATT	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/chain.sam" \
  --out "$tmpdir/chain" \
  --tmp-dir "$tmpdir/chain_buckets" \
  --threads 1 \
  --buckets 2 \
  --ham 1 \
  --ratio 0.60 \
  --emit-explain >/dev/null

grep -q 'AAAT	AAAA' "$tmpdir/chain.corrections.tsv"
grep -q 'AATT	AATT' "$tmpdir/chain.corrections.tsv"
if grep -q 'AATT	AAAA' "$tmpdir/chain.corrections.tsv"; then
  echo "chain correction should not exceed direct Hamming threshold" >&2
  exit 1
fi

cat > "$tmpdir/jitter.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:2000
j001	0	chr1	1000	60	10M100N10M	*	0	0	ACGTACGTAAACGTACGTAA	FFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
j002	0	chr1	1000	60	10M100N10M	*	0	0	ACGTACGTAAACGTACGTAA	FFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAA	UY:Z:IIII	GX:Z:GENE1
j003	0	chr1	1001	60	10M100N10M	*	0	0	ACGTACGTAAACGTACGTAA	FFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:AAAT	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/jitter.sam" \
  --out "$tmpdir/jitter" \
  --tmp-dir "$tmpdir/jitter_buckets" \
  --threads 1 \
  --buckets 2 \
  --ham 1 \
  --ratio 0.60 \
  --sj-jitter 10 \
  --emit-explain >/dev/null

grep -q 'AAAT	AAAA' "$tmpdir/jitter.corrections.tsv"

cat > "$tmpdir/structure.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:2000
s001	0	chr1	1000	60	10M100N10M	*	0	0	ACGTACGTAAACGTACGTAA	FFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE1
s002	0	chr1	1000	60	10M200N10M	*	0	0	ACGTACGTAAACGTACGTAA	FFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/structure.sam" \
  --out "$tmpdir/structure_default" \
  --tmp-dir "$tmpdir/structure_default_buckets" \
  --threads 1 \
  --buckets 2 \
  --emit-tsv >/dev/null

"$bin" \
  --bam "$tmpdir/structure.sam" \
  --out "$tmpdir/structure_baseline" \
  --tmp-dir "$tmpdir/structure_baseline_buckets" \
  --threads 1 \
  --buckets 2 \
  --no-structure \
  --emit-tsv >/dev/null

test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/structure_default.molecules.tsv")" = "2"
test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/structure_baseline.molecules.tsv")" = "1"

cat > "$tmpdir/rep.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
q_low	0	chr1	101	10	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:CCCC	UY:Z:IIII	GX:Z:GENE1
q_high	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELL1	UR:Z:CCCC	UY:Z:IIII	GX:Z:GENE1
SAM

"$bin" \
  --bam "$tmpdir/rep.sam" \
  --out "$tmpdir/rep" \
  --tmp-dir "$tmpdir/rep_buckets" \
  --threads 1 \
  --buckets 2 \
  --emit-tsv >/dev/null

grep -q '^q_low	.*	1	' "$tmpdir/rep.assignments.tsv"
grep -q '^q_high	.*	0	' "$tmpdir/rep.assignments.tsv"

cat > "$tmpdir/endbin.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:2000
e_pos_a	0	chr1	81	60	40M100N20M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE_END
e_pos_b	0	chr1	101	60	20M100N20M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE_END
e_neg_a	16	chr1	101	60	20M100N40M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE_END
e_neg_b	16	chr1	101	60	20M100N20M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF	CB:Z:CELL1	UR:Z:GGGG	UY:Z:IIII	GX:Z:GENE_END
SAM

"$bin" \
  --bam "$tmpdir/endbin.sam" \
  --out "$tmpdir/endbin_off" \
  --tmp-dir "$tmpdir/endbin_off_buckets" \
  --threads 1 \
  --buckets 2 \
  --emit-tsv >/dev/null

"$bin" \
  --bam "$tmpdir/endbin.sam" \
  --out "$tmpdir/endbin_on" \
  --tmp-dir "$tmpdir/endbin_on_buckets" \
  --threads 1 \
  --buckets 2 \
  --end-bin 50 \
  --emit-tsv >/dev/null

test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/endbin_off.molecules.tsv")" = "2"
test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/endbin_on.molecules.tsv")" = "4"

cat > "$tmpdir/multi_a.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
m_a	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELLX	UR:Z:TTTT	UY:Z:IIII	GX:Z:GENE_MULTI
SAM

cat > "$tmpdir/multi_b.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
m_b	0	chr1	101	60	10M	*	0	0	ACGTACGTAA	FFFFFFFFFF	CB:Z:CELLX	UR:Z:TTTT	UY:Z:IIII	GX:Z:GENE_MULTI
SAM

"$bin" \
  --bam "$tmpdir/multi_a.sam" \
  --bam "$tmpdir/multi_b.sam" \
  --out "$tmpdir/multi_shared" \
  --tmp-dir "$tmpdir/multi_shared_buckets" \
  --threads 1 \
  --buckets 2 \
  --emit-tsv >/dev/null

"$bin" \
  --bam "$tmpdir/multi_a.sam" \
  --bam "$tmpdir/multi_b.sam" \
  --out "$tmpdir/multi_isolated" \
  --tmp-dir "$tmpdir/multi_isolated_buckets" \
  --threads 1 \
  --buckets 2 \
  --isolate-inputs \
  --emit-tsv >/dev/null

test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/multi_shared.molecules.tsv")" = "1"
test "$(awk 'NR>1 {n++} END {print n+0}' "$tmpdir/multi_isolated.molecules.tsv")" = "2"
grep -q 'IN=input000000' "$tmpdir/multi_isolated.molecules.tsv"
grep -q 'IN=input000001' "$tmpdir/multi_isolated.molecules.tsv"
