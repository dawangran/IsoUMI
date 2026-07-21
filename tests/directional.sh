#!/bin/sh
set -eu

bin=${1:-src/isoumi}
tmpdir=$(mktemp -d "${TMPDIR:-/tmp}/isoumi-directional.XXXXXX")

cleanup() {
  rm -f "$tmpdir/input.sam"
  rm -f "$tmpdir/bad.stdout" "$tmpdir/bad.stderr"
  rm -f "$tmpdir/conflict.stdout" "$tmpdir/conflict.stderr"
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
