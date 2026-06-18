#pragma once
#include <htslib/sam.h>
char* build_group_key(bam1_t* b, const char* cell_tag, const char* gene_tag,
                      const char* source_tag, const char* input_scope_tag,
                      int no_gene, int no_structure, int locus_bin, int sj_jitter, int end_bin);
