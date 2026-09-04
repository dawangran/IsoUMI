# IsoUMI minimal WDL

Only `bams` is required. The minimal input file is:

```json
{
  "IsoUMIWorkflow.bams": [
    "/absolute/path/to/input.bam"
  ]
}
```

The workflow fixes the analysis choices used by this project:

- UMI-tools-style `directional` correction
- Hamming distance 1 per directional edge
- no gene-tag grouping
- no UMI-quality tag or quality-aware scoring
- chromosome, strand, splice-junction or locus structure retained
- corrected UMI written to `UB`
- duplicate state written to `DA` and SAM flag `0x400`
- molecule, assignment and correction reports always enabled
- representative reads always filtered, coordinate-sorted and indexed

Optional inputs and defaults are:

```text
output_prefix = "isoumi"
cell_tag = "CB"
raw_umi_tag = "UR"
isolate_inputs = false
end_bin = 0
threads = 4
buckets = 64
memory = "8 GB"
```

Disk allocation is intentionally not exposed by the workflow. Local runs use
the host filesystem, while cloud or HPC backends should configure task disk
space in their backend settings.

For multiple chunks of the same sample, list all BAMs and leave
`isolate_inputs=false`. For BAMs from different samples or libraries, set
`isolate_inputs=true` so reused cell barcodes cannot merge across inputs.

`end_bin=0` disables endpoint-aware grouping. Set a positive value such as 25
or 50 bp when transcript ends should also separate molecule groups.

Run with miniwdl:

```bash
miniwdl run wdl/IsoUMI.wdl -i wdl/example.inputs.json --dir workflow-output
```

Run with Cromwell:

```bash
java -jar cromwell.jar run wdl/IsoUMI.wdl --inputs wdl/example.inputs.json
```

The workflow uses `dawang02/isoumi:0.1.1` and produces both the full
duplicate-marked BAM and a sorted/indexed representative-read BAM.
