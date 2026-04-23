#pragma once
#include <stdint.h>
#include <htslib/sam.h>
int build_sj_or_locus(bam1_t* b, int locus_bin, int sj_jitter, int* is_spliced,
                      uint64_t* sj_hash_out, uint64_t* locus_hash_out);
