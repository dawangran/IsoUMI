# Directional UMI Correction Design

## Objective

Add an optional UMI-tools-compatible directional correction mode while preserving
IsoUMI's current direct count-ratio mode and isoform-aware grouping behavior.
The new mode must implement the directional graph rule, including transitive
network grouping, rather than merely replacing the existing ratio inequality.

## Approaches Considered

### 1. Replace the ratio inequality only

Change the direct raw-to-seed predicate from `raw / seed <= ratio` to
`seed >= 2 * raw - 1`.

This is small but is not UMI-tools directional correction because a UMI can only
reach the final seed directly. It misses valid decreasing-count chains whose
endpoints are farther apart than `--ham`.

### 2. Build directional networks inside each IsoUMI grouping key

Retain IsoUMI's existing cell, gene, strand, splice-junction/locus, input scope,
and optional transcript-end grouping. Within each resulting key, build directed
UMI edges using the UMI-tools rule and group each seed with all nodes reachable
through outgoing edges.

This is the selected approach. It combines IsoUMI's transcript-context isolation
with UMI-tools directional network semantics and needs no new dependency.

### 3. Invoke UMI-tools as an external process

Export each grouping key to UMI-tools and re-import the result. This would add a
Python/runtime dependency, complicate temporary I/O, and still require custom
translation of IsoUMI grouping keys and reports. It is rejected.

## Command-Line Interface

Add:

```text
--correction-method <ratio|directional>
```

The default is `ratio`, preserving all existing behavior and output for commands
that do not select the new mode.

- `ratio` uses the existing direct raw-to-seed algorithm and `--ratio`.
- `directional` uses Hamming distance plus `seed_count >= 2 * raw_count - 1`.
- `--ratio` remains accepted in directional mode for command compatibility but
  does not affect directional graph construction. Help and parameter
  documentation must state this explicitly.
- `--min-merge-confidence` is incompatible with `directional`. Parsing must fail
  with a clear message when a positive confidence threshold is combined with
  directional mode, because applying it would no longer reproduce the standard
  directional grouping rule.
- UMI quality may still be read for diagnostic report fields, but it must not
  affect directional edges, network membership, or seed ranking.

## Directional Algorithm

For the unique UMIs in one existing IsoUMI grouping key:

1. Sort candidates by descending count and then lexicographically. This gives a
   deterministic representative when counts tie.
2. For every equal-length UMI pair whose Hamming distance is at most `--ham`, add
   `A -> B` when `count(A) >= 2 * count(B) - 1`. Evaluate the reverse direction
   independently.
3. Visit candidate seeds in the deterministic order. If a seed has not already
   been assigned, traverse all outgoing edges with breadth-first search and
   assign every newly reached node to that seed.
4. Record a deterministic parent for every reached node. Neighbors are considered
   in the same count-descending, lexicographic order; the first discovered path
   owns the parent.
5. An unvisited UMI starts a new component and becomes its own seed.

This matches UMI-tools directional component membership while deliberately using
deterministic rather than random tie-breaking. The count expression must avoid
integer overflow by comparing in a wider integer type.

## Assignment and Reporting Model

Extend the internal assignment with:

- final representative index;
- immediate parent index;
- edge Hamming distance from the immediate parent;
- root Hamming distance from raw UMI to final representative;
- path length in graph edges;
- correction method.

The existing correction TSV columns retain their order and meaning. Append these
columns:

```text
method	parent_umi	parent_count	edge_hamming	path_length
```

Existing fields behave as follows:

- `corr_umi`, `seed_count`, and `seed_avgq` describe the final component root.
- `hamming` is the raw-to-final-root Hamming distance and may exceed `--ham` in
  directional mode.
- `confidence` remains `1.0` for self rows. For directional corrections it is a
  diagnostic heuristic computed for the selected immediate parent edge and is
  never used to accept or reject the edge.
- `reason` is `directional` for a directional correction and `self` for a root.
- Ratio-mode rows populate the appended fields consistently: the immediate
  parent is the root, `edge_hamming == hamming`, and path length is one; self rows
  use themselves as parent and path length zero.

No molecule or BAM schema changes are required. All reads assigned to a
directional component receive the root UMI as the corrected UMI.

## Error Handling and Compatibility

- Reject unknown correction methods during argument parsing.
- Reject positive `--min-merge-confidence` with directional mode regardless of
  option order.
- Continue rejecting unequal-length UMI pairs as potential neighbors; they stay
  in separate components unless linked through valid equal-length pairs, which
  cannot cross lengths.
- The default `ratio` mode must retain its direct-distance, ratio, quality, and
  confidence behavior.
- Existing scripts and report readers must continue to work because old columns
  are unchanged and new correction columns are appended.

## Tests

Add end-to-end shell assertions covering:

1. default mode remains `ratio` and keeps the current no-chain behavior;
2. `--correction-method directional` accepts the `3 -> 2` boundary that
   `--ratio 0.5` rejects;
3. a decreasing-count chain is one directional component even when the final
   raw-to-root Hamming distance exceeds `--ham`;
4. a count-invalid edge remains in a separate component;
5. singleton ties are deterministic;
6. report headers and parent/path fields describe both direct and transitive
   assignments;
7. unknown methods and directional plus positive confidence fail clearly;
8. all existing smoke, synthetic-truth, and BAM-tag tests remain green.

## Documentation

Update `README.md`, `docs/parameters.md`, the CLI help, and release notes to
describe both modes, the directional formula, transitive grouping, deterministic
tie-breaking, report-column semantics, and the fact that IsoUMI grouping context
still differs from a standalone UMI-tools run.

## Source Basis

The directional edge and connected-component semantics follow the UMI-tools
documentation and its `UMIClusterer._get_adj_list_directional` implementation:

- https://umi-tools.readthedocs.io/en/stable/the_methods.html
- https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/network.py

