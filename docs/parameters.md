# IsoUMI Parameter Guide

This guide explains how each IsoUMI parameter affects grouping, UMI correction,
deduplication, and reproducibility. The short help remains available with
`src/isoumi --help`.

## Inputs And Outputs

### `--bam <FILE>`

Input SAM/BAM file. This option is repeatable.

Use repeated `--bam` arguments when one run should combine multiple input files
before molecule correction:

```bash
src/isoumi --bam sample1.bam --bam sample2.bam --out merged
```

For multiple inputs, reference dictionaries must match by target name, order,
and length. The first available `M5` and `AS` values are preserved, and
conflicting values in any later input are rejected. Compatible metadata lines
are merged, while conflicting `@RG` IDs are rejected because they can change
sample/library semantics.

Program provenance is preserved by namespacing every `@PG` ID from inputs after
the first. IsoUMI rewrites `@PG ID`, `@PG PP`, and record-level `PG:Z:` values as
one mapping. For example, `minimap2` in the second input becomes
`minimap2.isoumi.input000001`. The output header also records the IsoUMI version
and invocation in a unique `@PG` entry. Every `PP` value must resolve to a
program ID declared in the same input header.

Use `--strip-pg` when program provenance is intentionally excluded to reduce
output size. In this mode IsoUMI removes every input `@PG` line, removes every
record-level `PG` auxiliary tag before bucket files are written, and does not
add its own `@PG` line. Other metadata, including `@RG`, `@SQ`, and `@CO`, is
unchanged. Removing the header and record tags together prevents dangling
`PG:Z:` references.

By default, repeated inputs are deduplicated together. This is appropriate for
technical splits or lanes from the same cell-barcode namespace. For independent
samples or libraries that may reuse cell barcodes, use `--isolate-inputs` or
include an existing source boundary with `--source-tag`.

### `--bam-list <FILE>`

Text file containing one SAM/BAM path per line. Blank lines and lines beginning
with `#` are ignored.

This is useful for large batches where repeated `--bam` arguments become hard
to read.

### `--out <PREFIX>`

Required output prefix. The corrected BAM is written as:

```text
<PREFIX>.dedup.bam
```

Optional reports use the same prefix:

```text
<PREFIX>.molecules.tsv
<PREFIX>.assignments.tsv
<PREFIX>.corrections.tsv
```

Final files are written through same-directory staging paths and atomically
renamed after a successful close. No final output may be the same path or inode
as an input BAM.


### `--strip-pg`

Remove all program-provenance metadata from temporary and final BAM files:

- all header `@PG` lines are omitted
- all record-level `PG` tags are deleted
- IsoUMI does not add its own version/command-line `@PG` line

This option does not change grouping, UMI correction, molecule selection,
duplicate marking, report contents, or other header records. Use it only when
the original command history is stored elsewhere and reduced BAM size is more
important than embedded provenance.

## Tag Parameters

All tag-name parameters must match `[A-Za-z][A-Za-z0-9]`. Tags assigned to
different semantic roles must be distinct. Invalid or conflicting tag
configurations fail before any output is created. `PG` is reserved for SAM
program provenance and cannot be selected by a tag-name option.

### `--cell-tag <TAG>` default `CB`

Cell barcode tag. Reads without this tag are passed through and do not
participate in UMI correction, molecule reports, or assignment reports.

Changing this is necessary for non-10x or custom preprocessing pipelines.

### `--umi-tag <TAG>` default `UR`

Raw UMI tag read from each input record. Reads without this tag are passed
through with duplicate flag `0`. Any pre-existing value in the configured
corrected-UMI or molecule-ID output tag is removed from those records.

UMI correction compares strings with Hamming distance, so UMIs must have equal
length to be compared.

### `--umi-qual-tag <TAG>` default `UY`

Optional UMI quality string. The string should have the same length as the UMI
and use Phred+33 characters.

When available, UMI quality contributes to ratio-mode seed ranking and
merge-confidence scoring. In directional mode it can contribute to the reported
diagnostic edge confidence, but never changes graph membership or seed ranking.
Missing or length-mismatched quality strings are ignored for that UMI.

### `--gene-tag <TAG>` default `GX`

Gene tag used in the grouping key unless `--no-gene` is set. Reads missing the
gene tag are grouped with `GX=NA`.

Use a gene tag that matches the annotation strategy used upstream.

### `--source-tag <TAG>`

Optional existing SAM tag to include in the grouping key. A common choice is
`RG` when read groups encode sample, library, or lane boundaries that should not
be deduplicated together.

Reads missing this tag are grouped with `SRC=NA`.

### `--input-scope-tag <TAG>` default `zi`

Temporary tag used internally by `--isolate-inputs`. The tag is added while
writing bucket BAMs, included in the grouping key as `IN=...`, and removed
before the final BAM is written.

The tag must not already exist in input records. If it does, choose a different
two-character tag.

### `--umi-out <TAG>` default `UB`

Output tag for the corrected UMI. For reads that are not corrected but have a
raw UMI, this tag is set to the raw UMI.

### `--dup-flag <TAG>` default `DA`

Integer duplicate flag written to each output record:

- `0`: representative read for a corrected molecule, or passthrough read
- `1`: duplicate read within a corrected molecule

Representative reads are selected by mapping quality, aligned reference span,
query length, and then stable input order.

IsoUMI also sets or clears the standard SAM duplicate bit (`0x400`) by default
so downstream BAM tools can recognize duplicates without reading `DA`.

Only primary alignments contribute to UMI counts and molecule representative
selection. Secondary (`0x100`) and supplementary (`0x800`) records inherit the
corrected UMI, duplicate status, and optional molecule tag from an unambiguous
matching primary record with the same QNAME, input scope, and read-end flag.
They are excluded from molecule and assignment TSV counts. An orphan
non-primary record passes through with `UB=UR` and duplicate status `0`.

Mapped secondary and supplementary records must carry the same cell tag as the
primary record. A mapped non-primary record without that tag is rejected because
it cannot be routed to the primary record's bucket safely.

### `--mol-tag <TAG>`

Optional output tag containing a molecule identifier:

```text
CB|grouping_key|corrected_UMI
```

This is useful when inspecting BAM records directly. The same identifier appears
in `*.assignments.tsv` when `--emit-tsv` is enabled.

Dynamic values percent-escape `%`, `|`, and `=`; a literal `NA` value is
also escaped so it cannot collide with the missing-value marker used in grouping
keys.

## Grouping Parameters

Grouping parameters decide which reads are allowed to compare UMIs. They do not
directly change UMI distance calculations.

### `--strand-mode <auto|alignment|ignore>` default `auto`

Controls whether BAM alignment direction (`FLAG 0x10`) separates molecule
groups. Alignment direction does not necessarily identify the RNA's transcription
strand. Reads expected to represent one RNA molecule may have opposite alignment
directions in some cDNA workflows before orientation normalization. Whether this
is expected depends on library preparation and upstream processing; it is not a
universal property of long-read RNA sequencing.

| Mode | Grouping behavior |
| --- | --- |
| `auto` | With gene grouping enabled and a nonempty, valid `Z` string in the selected `--gene-tag` (default `GX`), use `STR=.` and ignore alignment direction. Otherwise retain direction as `STR=+` or `STR=-`. |
| `alignment` | Always retain `FLAG 0x10` as `STR=+` or `STR=-`, reproducing the earlier grouping policy. |
| `ignore` | Always use `STR=.`, including missing gene tags and `--no-gene`. |

`auto` is a gene-tag heuristic, not protocol or transcript-strand inference. It
does not determine whether upstream processing has already oriented the reads.
Use `alignment` when direction should remain a molecule boundary, including
appropriately oriented inputs for which opposite alignments need to be kept
separate. The synthetic mixed-FLAG regression demonstrates grouping behavior;
it does not establish that a real sample contains mixed-direction molecules.

The decision in `auto` is made for each record. Missing, empty, or non-`Z` gene
tags retain alignment-direction separation. Every valid nonempty gene string is
trusted, including the literal `NA`; IsoUMI does not reinterpret placeholder or
multi-gene strings. Literal `NA` is escaped in grouping keys so it remains
distinct from a missing gene tag. Gene-tagged and untagged records still have
different gene keys and do not join merely because their other fields match.

This policy applies to exact and compatible structure grouping and to
`--no-structure`, with both ratio and directional UMI correction. Ignoring
direction allows otherwise matching reads in both directions to support one
molecule. Cell, gene, source/input scope, contig, selected structural context,
and UMI correction rules still define the group. It does not infer transcript
strand or reverse-complement UMI strings. Input sequence and the alignment
orientation flag are preserved; normal duplicate marking still applies.

**Migration:** `auto` is a default behavior change in the unreleased source.
Published IsoUMI 0.1.1 always separates alignment directions. To reproduce that
grouping on current source, add `--strand-mode alignment` and keep the other
analysis settings unchanged. Gene-tagged reads now use different grouping keys
and molecule IDs under `auto`, even when all reads align in one direction and
the molecule count is unchanged. With `--end-bin`, ignored direction also uses
genomic `EL`/`ER` terms instead of the legacy `E5`/`E3` terms.

`auto` relies on upstream annotations to distinguish genes and antisense
features that should remain separate. It cannot distinguish same-gene antisense
molecules assigned the same gene string; shared UMIs can still collide.
`ignore` broadens grouping without requiring gene annotations, so choose it
only when that interpretation fits the input. This policy does not establish
biological molecule identity by itself.

### `--structure-mode <exact|compatible>` default `exact`

`exact` uses the existing SJ grid and locus-bin keys. `compatible` is an
experimental alternative for variable alignment coordinates and truncated
junction chains. It retains cell, gene, source, input, contig and the direction
boundaries selected by `--strand-mode`. Within each boundary it uses
equal-length raw UMI Hamming neighborhoods to limit the search, compares
structures using the original coordinates, then reruns the selected UMI
correction method using counts in the final structure groups. A candidate
neighborhood is not itself a molecule.

Observations are ordered by junction count, reference span and coordinates,
with deterministic ties. Partial reads can attach to compatible longer
structures. Fixed structural witnesses and checks against assigned members
prevent explicit splice/retained-intron conflicts from joining the same group;
lack of overlap between two partial reads is not itself a conflict. An anchor
still requires positive overlap with a candidate read. Structural compatibility
does not connect conflicting anchors transitively.

Anchor support is frozen before ambiguous observations are assigned. If all
otherwise compatible anchors conflict with already assigned members, the
observation remains in a separate unresolved group. Low support does not justify
forcing a conflicting observation into a dominant group. Fixed witnesses make
this a conservative heuristic and can retain extra groups. The assignment
report exposes ambiguities; it does not establish a uniquely identified source
molecule.

This mode cannot be combined with `--no-structure` or `--end-bin`.
`--sj-jitter` and `--locus-bin` apply only to exact mode.

### `--sj-tolerance <INT>` default `10`

Maximum absolute difference in base pairs when matching each boundary of an
observed junction to a compatible-mode anchor junction. Zero requires exact
boundary coordinates. This is a direct coordinate tolerance, distinct from
`--sj-jitter` grid rounding. Increasing the tolerance can accommodate alignment
variation but also removes the ability to distinguish genuine nearby splice
sites. Truncation compatibility must still satisfy the junction-chain and
exonic-coverage checks.

### `--min-structure-support <INT>` default `3`

Primary-read support threshold used to prefer anchors when a read is compatible
with more than one structure. Among equally preferred anchors, raw-UMI support
and then total anchor support determine the choice. A uniquely compatible read
can join an anchor below this threshold. Conflicting low-support structures
remain separate. This heuristic is not a confidence level or proof of molecular
identity; it applies only to compatible mode.

### `--no-gene`

Ignore the gene tag during grouping. The grouping key uses `GX=NA` for all
reads. With the default `--strand-mode auto`, this retains alignment-direction
separation. Add `--strand-mode ignore` to omit that boundary as well.

Use this only when gene tags are absent, unreliable, or intentionally excluded
from an analysis. Removing the gene boundary can increase UMI comparisons;
under `auto`, it also restores the alignment-direction boundary.

### `--no-structure`

Baseline mode for application-note comparisons. It ignores splice-junction,
locus, and transcript-end structure in the grouping key.

With this option, grouping still uses:

```text
cell barcode + gene/NA + optional source/input scope + reference target + strand-mode policy
```

This mode is intentionally less isoform-aware and is useful for showing how much
splice-junction and locus context prevent over-collapse.

### `--isolate-inputs`

Treat each `--bam` argument or `--bam-list` entry as a separate source during
grouping. IsoUMI writes an internal input-scope tag into temporary bucket BAMs,
adds it to the molecule key, and removes it before writing the final BAM.

Use this when separate samples or libraries may reuse the same cell barcode
values. Do not use it for technical splits that should be deduplicated together
across files.

### `--locus-bin <INT>` default `1000`

Bin size in base pairs for non-spliced reads, defined as reads without a CIGAR
`N` operation.

For non-spliced reads, IsoUMI bins the reference start and end positions:

```text
start_bin = floor(start / locus_bin) * locus_bin
end_bin   = floor(end   / locus_bin) * locus_bin
```

The pair is hashed into an `LX=...` grouping component.

Smaller values are more precise and reduce over-grouping. Larger values tolerate
more alignment-end variation but can merge nearby contexts.

Typical choices:

- `1000`: default broad locus grouping for noisy long-read alignments
- `100` or `250`: stricter locus separation
- `2000`: more permissive grouping for sparse data

### `--sj-jitter <INT>` default `10`

Rounding window for splice-junction boundaries before hashing. It applies only
to spliced reads with CIGAR `N` operations.

Each junction start/end is rounded to the nearest multiple of `sj-jitter`.
For example, with `--sj-jitter 10`:

```text
1004 -> 1000
1005 -> 1010
1009 -> 1010
```

Lower values are stricter. Higher values tolerate small alignment boundary
differences but can merge nearby alternative splice sites.

Typical choices:

- `0` or `1`: exact junction boundaries
- `10`: default tolerance for small long-read alignment jitter
- `20` or `50`: more permissive; inspect results carefully

### `--end-bin <INT>`

Adds alignment-end bins to the grouping key in exact structure mode. The end
labels follow the direction policy chosen by `--strand-mode`.

When alignment direction is ignored (`STR=.`), the bins use genomic left and
right coordinates, independently of `FLAG 0x10`:

```text
EL = read_start
ER = read_end
```

When alignment direction is retained, forward-aligned reads use:

```text
E5 = read_start
E3 = read_end
```

Reverse-aligned reads use:

```text
E5 = read_end
E3 = read_start
```

Both ends are binned with `floor(position / end_bin) * end_bin`; `read_end` is
exclusive. The `E5`/`E3` labels follow BAM alignment direction and do not infer
RNA transcription strand. `--strand-mode alignment` preserves the legacy end
keys. The genomic `EL`/`ER` keys prevent opposite alignment directions from
splitting an otherwise matching group when direction is ignored.

Use this when transcript-end differences are biologically important. Smaller
values are more isoform-specific but more sensitive to truncation, soft
clipping, and alignment variation.

Typical choices:

- disabled: default, groups by junction or locus structure only
- `25` or `50`: transcript-end-aware grouping
- `100` or `200`: softer end-aware grouping for noisy data

## UMI Correction Parameters

UMI correction is performed inside each grouping key. Select direct ratio
assignment or directional graph grouping with `--correction-method`.

### `--correction-method <ratio|directional>` default `ratio`

`ratio` preserves the original IsoUMI behavior. Candidate seeds are processed
from strongest to weakest, and each raw UMI must directly satisfy the Hamming,
count-ratio, and optional confidence filters against its final seed.

`directional` builds UMI-tools-style directed networks. For two UMIs within the
Hamming threshold, an edge from a potential parent to a child is allowed when:

```text
seed_count >= 2 * raw_count - 1
```

Each unassigned UMI is visited in descending count and lexicographic order and
becomes the root of all nodes reachable through outgoing edges. This permits
transitive corrections through decreasing-count chains. The grouping context is
still IsoUMI's cell, gene, direction policy, and splice-junction or locus
context, so results are not identical to a standalone UMI-tools run configured
with different grouping.

### `--ham <INT>` default `1`

Maximum Hamming distance used for a direct ratio correction or for each edge in
a directional network.

In ratio mode, `--ham 1` means the final reported `raw_umi -> corr_umi`
correction differs by at most one position. In directional mode each edge must
differ by at most one position, but a multi-edge path can produce a final
raw-to-root `hamming` value greater than one.

Use `0` to disable sequence-error correction while still marking exact-UMI
duplicates.

### `--ratio <FLOAT>` default `0.10`, range `0..1`

Maximum allowed support ratio for `--correction-method ratio`:

```text
raw_count / seed_count <= ratio
```

Lower values are more conservative. Higher values correct more aggressively and
can merge true low-abundance molecules if the grouping context is broad.

Typical choices:

- `0.05`: conservative
- `0.10`: default
- `0.50` or `0.60`: useful for small tests or highly error-prone UMIs
- `1.00`: permits equal-count correction according to seed ranking

Directional mode accepts this option for command compatibility but ignores it.

### `--min-merge-confidence <FLOAT>` default `0.00`, range `0..1`

Optional heuristic confidence floor for ratio-mode corrections. The score
combines count ratio, UMI distance, and optional quality evidence. It is not a
posterior probability.

Set this above zero when you want a stricter correction policy:

```bash
--min-merge-confidence 0.45
--min-merge-confidence 0.60
```

A positive value is incompatible with `--correction-method directional` because
filtering graph edges by this score would no longer implement the standard
directional count/network rule. Directional reports still include diagnostic
confidence values when `--emit-explain` is enabled.

### `--no-quality-aware`

Disable quality-aware ratio-mode seed ranking and quality contribution to merge
confidence. Directional graph membership is quality-independent with or without
this option.

Use this when UMI quality tags are absent, unreliable, or not comparable across
inputs.

## Performance And Temporary Files

### `--threads <INT>` default `4`

Number of worker threads for per-bucket processing. Threading happens across
cell-barcode buckets; the same setting is also used for parallel input-BAM
decoding during sharding and parallel compression of the final BAM.

If IsoUMI is built without OpenMP support, this option is accepted but
per-bucket processing runs serially. The program prints a runtime notice in
that case.

### `--buckets <INT>` default `64`

Number of cell-barcode buckets. Increasing this can lower peak memory per bucket
on large datasets, but creates more temporary files.

Within each bucket, IsoUMI stream-aggregates repeated grouping-key/UMI pairs and
hash-deduplicates mapped non-primary read identities. Peak statistical working
memory therefore follows the number of unique pairs and unique non-primary
identities in simultaneously active buckets, not the raw record count. A single
extremely high-complexity cell can still dominate one bucket.

Typical choices:

- `16` or `32`: small datasets
- `64`: default
- `128` or `256`: larger datasets or high-depth cells

For hundred-GB-scale BAMs, a practical starting point is 16 CPU threads, 64 GB
RAM, and 128-256 buckets. IsoUMI does not use a GPU. Increase RAM or bucket count
after measuring a representative sample because unique molecule complexity and
cell imbalance, rather than compressed BAM size alone, determine peak memory.

### `--tmp-dir <DIR>` default `<out>.isoumi.tmp.<pid>`

Directory for bucket BAMs and temporary per-bucket reports. The directory must
not already exist; this prevents accidental reuse of stale bucket files.

Intermediate BAMs use fast level-1 compression because they are read only by
later pipeline phases and recompressed into the final BAM. For hundred-GB-scale
inputs, use a local SSD/NVMe filesystem with enough free space; a network
filesystem can make sharding and concatenation I/O-bound.

By default, temporary BAMs are reclaimed incrementally. A source bucket is
removed only after its deduplicated replacement closes successfully. During
final concatenation, each deduplicated bucket and TSV fragment is reclaimed
after it has been copied successfully into a staging output. The final BAM and
reports are atomically renamed into place only after a successful close, so an
existing output cannot be left partially overwritten. Consequently, normal peak
scratch usage is approximately one intermediate BAM generation plus the growing
staged final output, report data, and in-flight bucket overhead. Compression
ratios and uneven bucket sizes still matter. The output directory also needs
space for the staging file even when `--tmp-dir` is on another filesystem.

### `--keep-tmp`

Keep temporary bucket files after the run. Use this for debugging bucket-level
failures or inspecting intermediate BAMs. It disables incremental BAM deletion,
so both source and deduplicated bucket generations remain on disk and require
substantially more temporary space.

## Report Parameters

### `--emit-tsv`

Write molecule and assignment reports:

```text
<out>.molecules.tsv
<out>.assignments.tsv
```

Use this for benchmarking, molecule count comparisons, and debugging duplicate
labels.

Compatible mode adds a `structure_status` column to the assignment report:

| Value | Meaning |
| --- | --- |
| `compatible` | A unique compatible anchor whose unambiguous support reaches `--min-structure-support` |
| `unsupported` | An anchor below the support threshold, or a separate unresolved group created after conflict checks; the observation is retained |
| `ambiguous` | More than one compatible anchor; assignment uses frozen support and deterministic tie breaking |

These labels describe structural assignment, not calibrated probabilities.
Molecule counts and duplicate flags still come from local UMI correction.

### `--emit-explain`

Write correction explanations:

```text
<out>.corrections.tsv
```

This is the most useful report for parameter tuning. It shows raw UMI,
corrected UMI, counts, Hamming distance, quality summaries, confidence, reason,
bucket, method, immediate parent, edge distance, and graph path length.

The correction report uses `corr_umi`, `seed_count`, and `seed_avgq` for the
final component root before lower-support UMIs are merged into it.
`parent_umi`, `parent_count`, and `edge_hamming` describe the immediate edge;
`path_length` distinguishes self, direct, and transitive assignments. Final
molecule read counts are in `<out>.molecules.tsv`.

### `--no-bam-dup-flag`

Leave the standard SAM duplicate bit unchanged. The custom duplicate tag
selected by `--dup-flag` is still written.

## Recommended Profiles

### Default exploratory run

```bash
src/isoumi --bam input.bam --out sample --emit-tsv --emit-explain
```

### More conservative correction

```bash
src/isoumi \
  --bam input.bam \
  --out conservative \
  --ham 1 \
  --ratio 0.05 \
  --min-merge-confidence 0.60 \
  --emit-explain
```

### Directional correction

```bash
src/isoumi \
  --bam input.bam \
  --out directional \
  --correction-method directional \
  --ham 1 \
  --emit-tsv \
  --emit-explain
```

### Alignment-end-aware grouping

```bash
src/isoumi \
  --bam input.bam \
  --out end_aware \
  --end-bin 50 \
  --emit-tsv \
  --emit-explain
```

### Application-note baseline

```bash
src/isoumi \
  --bam input.bam \
  --out baseline_no_structure \
  --no-structure \
  --emit-tsv \
  --emit-explain
```
