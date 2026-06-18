#pragma once
#include <stddef.h>
#include "vector.h"

typedef struct {
  strvec_t bam_list;
  char *out_prefix;
  char *tmp_dir;
  int buckets, threads;
  char *cell_tag, *umi_tag, *umi_qual_tag, *gene_tag, *source_tag;
  char *umi_out, *dup_flag, *mol_tag, *input_scope_tag;
  int  no_gene, no_structure, ham, locus_bin, sj_jitter, end_bin, emit_tsv, emit_explain, quality_aware;
  int  keep_tmp, isolate_inputs, set_bam_dup_flag;
  double ratio, min_merge_confidence;
} cli_opts_t;

int parse_args(int argc, char **argv, cli_opts_t *o);
void free_opts(cli_opts_t *o);
