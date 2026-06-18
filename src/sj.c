#include "sj.h"
#include "key.h"
#include <stdio.h>
#include <stdint.h>

static inline uint64_t hash_u32pair(uint64_t h, uint32_t a, uint32_t b){
  const unsigned char* p = (const unsigned char*)&a;
  for (int i=0;i<4;++i){ h ^= p[i]; h *= 1099511628211ull; }
  p = (const unsigned char*)&b;
  for (int i=0;i<4;++i){ h ^= p[i]; h *= 1099511628211ull; }
  return h;
}
static inline int round_by_jitter(int x, int j){
  if (j<=1) return x;
  if (x>=0) return ((x + j/2) / j) * j;
  return ((x - j/2) / j) * j;
}

int build_sj_or_locus(bam1_t* b, int locus_bin, int sj_jitter, int* is_spliced,
                      uint64_t* sj_hash_out, uint64_t* locus_hash_out) {
  *is_spliced = 0; *sj_hash_out = 0; *locus_hash_out = 0;

  if (b->core.n_cigar == 0 || (b->core.flag & BAM_FUNMAP)) {
    int32_t pos = b->core.pos;
    int32_t end_pos = bam_endpos(b);
    int32_t sbin = (pos / locus_bin) * locus_bin;
    int32_t ebin = (end_pos / locus_bin) * locus_bin;
    uint64_t h = 1469598103934665603ull;
    h = hash_u32pair(h, (uint32_t)sbin, (uint32_t)ebin);
    *locus_hash_out = h;
    return 0;
  }

  uint32_t* cigar = bam_get_cigar(b);
  int32_t refpos = b->core.pos;

  uint64_t h = 1469598103934665603ull;
  int has_N = 0;

  for (int i=0; i < (int)b->core.n_cigar; ++i) {
    int op  = bam_cigar_op(cigar[i]);
    int len = bam_cigar_oplen(cigar[i]);
    if (op == BAM_CREF_SKIP) {
      has_N = 1;
      int32_t start = refpos;
      int32_t end   = refpos + len;
      if (sj_jitter>1){ start = round_by_jitter(start, sj_jitter); end = round_by_jitter(end, sj_jitter); }
      h = hash_u32pair(h, (uint32_t)start, (uint32_t)end);
      refpos += len;
    } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
      refpos += len;
    }
  }

  if (has_N) {
    *is_spliced = 1;
    *sj_hash_out = h;
  } else {
    int32_t pos = b->core.pos;
    int32_t end_pos = bam_endpos(b);
    int32_t sbin = (pos / locus_bin) * locus_bin;
    int32_t ebin = (end_pos / locus_bin) * locus_bin;
    uint64_t hx = 1469598103934665603ull;
    hx = hash_u32pair(hx, (uint32_t)sbin, (uint32_t)ebin);
    *locus_hash_out = hx;
  }
  return 0;
}
