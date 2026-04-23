#pragma once
#include <htslib/sam.h>
char* build_group_key(bam1_t* b, const char* cell_tag, const char* gene_tag,
                      int no_gene, int locus_bin, int sj_jitter, int end_bin);
