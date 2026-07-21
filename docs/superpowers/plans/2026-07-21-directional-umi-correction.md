# Directional UMI Correction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an optional UMI-tools-compatible directional graph correction mode without changing IsoUMI's default direct-ratio behavior.

**Architecture:** Parse a correction-method enum in the existing CLI, then dispatch each grouping key to either the current direct seed assignment or a new sparse-adjacency directional breadth-first traversal. Extend assignment metadata so the correction report can distinguish the final root from the immediate graph parent while leaving BAM and molecule output schemas unchanged.

**Tech Stack:** C11, htslib, POSIX shell integration tests, GNU Make.

## Global Constraints

- Default correction method remains `ratio`.
- Directional edges use `count(A) >= 2 * count(B) - 1` and Hamming distance `<= --ham`.
- Directional mode ignores `--ratio` and rejects positive `--min-merge-confidence`.
- UMI quality must not affect directional graph membership or representative ranking.
- Existing correction TSV columns retain their order; new columns are appended.
- No external UMI-tools or Python runtime dependency may be added.
- Directional representative ties are deterministic: count descending, then UMI lexicographic order.

---

### Task 1: Correction-method CLI contract

**Files:**
- Create: `tests/directional.sh`
- Modify: `src/Makefile`
- Modify: `src/cli.h`
- Modify: `src/cli.c`

**Interfaces:**
- Consumes: existing `parse_args(int argc, char **argv, cli_opts_t *o)`.
- Produces: `correction_method_t`, `cli_opts_t.correction_method`, and `--correction-method <ratio|directional>`.

- [ ] **Step 1: Write the failing CLI test**

Create `tests/directional.sh` with a temporary one-read SAM and assertions that:

```sh
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
```

Add `sh ../tests/directional.sh ./$(BIN)` to `src/Makefile` after the existing tests.

- [ ] **Step 2: Run the CLI test to verify RED**

Run:

```sh
make -C src test
```

Expected: existing tests pass, then `tests/directional.sh` fails because the help text lacks `--correction-method` and the parser rejects it as an unknown option.

- [ ] **Step 3: Add the correction-method enum and parser**

Add to `src/cli.h`:

```c
typedef enum {
  CORRECTION_RATIO = 0,
  CORRECTION_DIRECTIONAL = 1
} correction_method_t;
```

Add `correction_method_t correction_method;` to `cli_opts_t`.

In `src/cli.c`, default the field to `CORRECTION_RATIO`, add long option id `29`, and parse exactly `ratio` and `directional`:

```c
static int parse_correction_method(const char *value, correction_method_t *out){
  if (strcmp(value, "ratio") == 0){
    *out = CORRECTION_RATIO;
    return 0;
  }
  if (strcmp(value, "directional") == 0){
    *out = CORRECTION_DIRECTIONAL;
    return 0;
  }
  fprintf(stderr, "unknown --correction-method: %s\n", value);
  return -1;
}
```

After option parsing, reject the incompatible positive confidence setting:

```c
if (o->correction_method == CORRECTION_DIRECTIONAL &&
    o->min_merge_confidence > 0.0){
  fprintf(stderr,
          "--min-merge-confidence is incompatible with --correction-method directional\n");
  return -1;
}
```

Add the exact help entry:

```text
  --correction-method <ratio|directional>
                         UMI correction algorithm (default: ratio)
```

- [ ] **Step 4: Run the CLI test to verify GREEN**

Run:

```sh
make -C src test
```

Expected: all existing tests and the new CLI assertions pass.

- [ ] **Step 5: Commit the CLI contract**

```sh
git add tests/directional.sh src/Makefile src/cli.h src/cli.c
git commit -m "feat: add directional correction CLI mode"
```

---

### Task 2: Directional graph assignment and audit metadata

**Files:**
- Modify: `tests/directional.sh`
- Modify: `src/pipeline.c`

**Interfaces:**
- Consumes: `cli_opts_t.correction_method`, sorted `umi_stat_t` arrays, and `--ham`.
- Produces: ratio or directional `umi_assignment_t` values containing root, parent, edge distance, root distance, path length, and diagnostic confidence.

- [ ] **Step 1: Write failing directional behavior tests**

Extend `tests/directional.sh` with one SAM containing four independent cell groups:

```text
CELL_CHAIN:    AAAA=5, AAAT=3, AATT=2, ATTT=2
CELL_BOUNDARY: CCCC=3, CCCT=2
CELL_TIE:      GGGA=1, GGGG=1
CELL_DIRECT:   TTTT=10, TTTA=1
```

Run ratio and directional outputs with `--ham 1 --emit-tsv --emit-explain`.
Assert the default ratio run keeps `AATT` separate and the directional run has:

```sh
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
```

Also assert the header is exactly the old eleven columns followed by:

```text
method	parent_umi	parent_count	edge_hamming	path_length
```

For `CELL_DIRECT`, assert the ratio row reports method `ratio`, parent `TTTT`,
parent count `10`, edge Hamming `1`, and path length `1`.

- [ ] **Step 2: Run the directional behavior test to verify RED**

Run:

```sh
sh tests/directional.sh src/isoumi
```

Expected: FAIL because directional currently dispatches to the ratio algorithm and the appended report fields do not exist.

- [ ] **Step 3: Extend assignment and report structures**

In `src/pipeline.c`, extend `umi_assignment_t`:

```c
typedef struct {
  int rep_idx;
  int parent_idx;
  int hamming;
  int edge_hamming;
  int path_length;
  double corr_mismatch_q;
  double raw_mismatch_q;
  double confidence;
  int mismatch_q_n;
} umi_assignment_t;
```

Extend `map_item_t` with `char *parent_umi`, `int parent_count`,
`int edge_hamming`, `int path_length`, and `correction_method_t method`. Update
`free_map_items` to free `parent_umi`.

Initialize ratio self rows with parent=self and path zero. For a direct ratio
assignment, set parent=root, edge Hamming equal to root Hamming, and path one.

Append the new fields in `emit_correction_row` and change `merge_reason` so a
non-self directional row reports `directional`.

- [ ] **Step 4: Implement sparse directional graph traversal**

Add a local integer-vector adjacency type:

```c
typedef struct {
  int *data;
  int n;
  int cap;
} int_vec_t;
```

Provide `int_vec_push` and `free_directional_graph` helpers. Build pairwise edges
only for equal-length UMIs with Hamming distance at most `o->ham`; compare counts
using `int64_t`:

```c
static int directional_count_allows(const umi_stat_t *from,
                                    const umi_stat_t *to){
  return (int64_t)from->count >= (int64_t)2 * (int64_t)to->count - 1;
}
```

Implement:

```c
static int build_directional_assignments(const cli_opts_t *o,
                                         const umi_stat_t *stats,
                                         int n,
                                         umi_assignment_t **out);
```

It allocates the sparse adjacency list and an `n`-element queue, visits seeds in
the already sorted array order, and performs breadth-first traversal. On first
discovery of child `v` from parent `u`, populate:

```c
assignments[v].rep_idx = seed;
assignments[v].parent_idx = u;
assignments[v].edge_hamming = edge_ham;
assignments[v].path_length = assignments[u].path_length + 1;
assignments[v].hamming = root_ham;
assignments[v].confidence = compute_merge_confidence(
    o, &stats[u], &stats[v], edge_ham, parent_mq, raw_mq, mismatch_q_n);
```

Compute `root_ham` with `compute_umi_distance_metrics(..., -1, ...)` so a valid
chain can report a final-root distance larger than `--ham`. Dispatch through:

```c
static int build_umi_assignments(const cli_opts_t *o,
                                 const umi_stat_t *stats,
                                 int n,
                                 umi_assignment_t **out){
  if (o->correction_method == CORRECTION_DIRECTIONAL){
    return build_directional_assignments(o, stats, n, out);
  }
  return build_ratio_assignments(o, stats, n, out);
}
```

- [ ] **Step 5: Populate report metadata and verify GREEN**

When translating assignments to `map_item_t`, duplicate the parent UMI and copy
the parent count, edge distance, path length, and method. Set
`quality_supported` only in ratio mode.

Change the correction header to:

```text
key	raw_umi	corr_umi	raw_count	seed_count	hamming	raw_avgq	seed_avgq	confidence	reason	bucket	method	parent_umi	parent_count	edge_hamming	path_length
```

Run:

```sh
make -C src test
```

Expected: all directional boundary, chain, invalid-edge, tie, and report tests pass along with the existing suite.

- [ ] **Step 6: Commit the directional engine**

```sh
git add tests/directional.sh src/pipeline.c
git commit -m "feat: implement directional UMI networks"
```

---

### Task 3: User documentation and release verification

**Files:**
- Modify: `README.md`
- Modify: `docs/parameters.md`
- Modify: `docs/application_note.md`
- Modify: `CHANGELOG.md`

**Interfaces:**
- Consumes: the implemented CLI and report schema.
- Produces: accurate user-facing method selection, semantics, limitations, and report documentation.

- [ ] **Step 1: Write the documentation acceptance checks**

Extend `tests/directional.sh` to assert CLI help names both methods. Extend
`scripts/check_release.sh` with literal checks for `--correction-method`,
`seed_count >= 2 * raw_count - 1`, and `parent_umi` in the primary documentation.

- [ ] **Step 2: Run release checks to verify RED**

Run:

```sh
sh scripts/check_release.sh
```

Expected: FAIL because the primary documentation does not yet describe the new mode and report fields.

- [ ] **Step 3: Update method and report documentation**

Document in `README.md` and `docs/parameters.md`:

```text
--correction-method ratio
--correction-method directional
```

State that ratio is the backward-compatible default, directional ignores
`--ratio`, directional rejects positive `--min-merge-confidence`, graph paths may
produce raw-to-root Hamming distances above `--ham`, and IsoUMI grouping remains
isoform-aware rather than identical to standalone UMI-tools grouping.

Document appended correction columns and distinguish `corr_umi` from
`parent_umi`. Update `docs/application_note.md` so its implementation and
limitations describe both modes. Add an `Unreleased` entry to `CHANGELOG.md`.

- [ ] **Step 4: Run documentation and full verification**

Run:

```sh
make clean
make
make test
git diff --check
```

Expected: clean rebuild succeeds with no compiler warnings, all shell tests and release checks pass, and `git diff --check` prints no output.

- [ ] **Step 5: Commit documentation and checks**

```sh
git add README.md docs/parameters.md docs/application_note.md CHANGELOG.md \
  scripts/check_release.sh tests/directional.sh
git commit -m "docs: document directional UMI correction"
```

