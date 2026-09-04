version 1.0

task RunIsoUMI {
  input {
    Array[File]+ bams
    String output_prefix = "isoumi"
    String cell_tag = "CB"
    String raw_umi_tag = "UR"
    Boolean isolate_inputs = false
    Int end_bin = 0
    Int threads = 4
    Int buckets = 64
    String memory = "32 GB"
  }

  command <<<
    set -eu

    out_prefix="~{output_prefix}"
    set -- \
      --bam '~{sep="' --bam '" bams}' \
      --out "$out_prefix" \
      --cell-tag "~{cell_tag}" \
      --umi-tag "~{raw_umi_tag}" \
      --umi-out UB \
      --dup-flag DA \
      --mol-tag MI \
      --correction-method directional \
      --ham 1 \
      --no-gene \
      --no-quality-aware \
      --threads "~{threads}" \
      --buckets "~{buckets}" \
      --emit-tsv \
      --emit-explain

    if [ "~{isolate_inputs}" = "true" ]; then set -- "$@" --isolate-inputs; fi
    if [ "~{end_bin}" -gt 0 ]; then set -- "$@" --end-bin "~{end_bin}"; fi

    if ! isoumi "$@"; then
      echo "IsoUMI failed; representative-read outputs were not attempted" >&2
      exit 1
    fi

    unsorted_bam="$out_prefix.representatives.unsorted.bam"
    samtools view -@ "~{threads}" -b -F 3332 -d MI \
      -o "$unsorted_bam" "$out_prefix.dedup.bam"
    samtools sort -@ "~{threads}" \
      -o "$out_prefix.representatives.bam" "$unsorted_bam"
    rm -f "$unsorted_bam"
    samtools index -@ "~{threads}" "$out_prefix.representatives.bam"
  >>>

  output {
    File duplicate_marked_bam = output_prefix + ".dedup.bam"
    File representatives_bam = output_prefix + ".representatives.bam"
    File representatives_bai = output_prefix + ".representatives.bam.bai"
    File molecules_tsv = output_prefix + ".molecules.tsv"
    File assignments_tsv = output_prefix + ".assignments.tsv"
    File corrections_tsv = output_prefix + ".corrections.tsv"
  }

  runtime {
    docker: "dawang02/isoumi:0.1.1"
    cpu: threads
    memory: memory
  }

  meta {
    description: "Minimal annotation-independent IsoUMI workflow using directional UMI correction with Hamming distance 1."
  }
}

workflow IsoUMIWorkflow {
  input {
    Array[File]+ bams
    String output_prefix = "isoumi"
    String cell_tag = "CB"
    String raw_umi_tag = "UR"
    Boolean isolate_inputs = false
    Int end_bin = 0
    Int threads = 4
    Int buckets = 64
    String memory = "32 GB"
  }

  call RunIsoUMI {
    input:
      bams = bams,
      output_prefix = output_prefix,
      cell_tag = cell_tag,
      raw_umi_tag = raw_umi_tag,
      isolate_inputs = isolate_inputs,
      end_bin = end_bin,
      threads = threads,
      buckets = buckets,
      memory = memory
  }

  output {
    File duplicate_marked_bam = RunIsoUMI.duplicate_marked_bam
    File representatives_bam = RunIsoUMI.representatives_bam
    File representatives_bai = RunIsoUMI.representatives_bai
    File molecules_tsv = RunIsoUMI.molecules_tsv
    File assignments_tsv = RunIsoUMI.assignments_tsv
    File corrections_tsv = RunIsoUMI.corrections_tsv
  }
}
