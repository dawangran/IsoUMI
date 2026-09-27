#pragma once
#include <stddef.h>
#include "vector.h"

#define ISOUMI_VERSION_NUMBER "0.1.1"

typedef enum {
  CORRECTION_RATIO = 0,
  CORRECTION_DIRECTIONAL = 1
} correction_method_t;

typedef enum { STRUCTURE_EXACT = 0, STRUCTURE_COMPATIBLE = 1 } structure_mode_t;

typedef enum {
  STRAND_AUTO = 0,
  STRAND_ALIGNMENT = 1,
  STRAND_IGNORE = 2
} strand_mode_t;

typedef struct {
  strvec_t bam_list;
  char *out_prefix, *command_line;
  char *tmp_dir;
  int buckets, threads;
  char *cell_tag, *umi_tag, *umi_qual_tag, *gene_tag, *source_tag;
  char *umi_out, *dup_flag, *mol_tag, *input_scope_tag;
  int  no_gene, no_structure, ham, locus_bin, sj_jitter, end_bin, emit_tsv, emit_explain, quality_aware;
  int  keep_tmp, isolate_inputs, set_bam_dup_flag, strip_pg;
  correction_method_t correction_method;
  structure_mode_t structure_mode;
  strand_mode_t strand_mode;
  int sj_tolerance, min_structure_support;
  double ratio, min_merge_confidence;
} cli_opts_t;

int parse_args(int argc, char **argv, cli_opts_t *o);
void free_opts(cli_opts_t *o);
