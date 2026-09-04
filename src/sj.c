#include "sj.h"
#include "key.h"
#include <stdio.h>
#include <stdint.h>

static inline uint64_t hash_u64pair(uint64_t h, uint64_t a, uint64_t b){
  for (unsigned int shift=0; shift<64; shift+=8){
    h ^= (unsigned char)(a >> shift);
    h *= 1099511628211ull;
  }
  for (unsigned int shift=0; shift<64; shift+=8){
    h ^= (unsigned char)(b >> shift);
    h *= 1099511628211ull;
  }
  return h;
}
static inline hts_pos_t round_by_jitter(hts_pos_t x, int j){
  if (j<=1) return x;
  if (x>=0) return ((x + j/2) / j) * j;
  return ((x - j/2) / j) * j;
}

int build_sj_or_locus(bam1_t* b, int locus_bin, int sj_jitter, int* is_spliced,
                      uint64_t* sj_hash_out, uint64_t* locus_hash_out) {
  *is_spliced = 0; *sj_hash_out = 0; *locus_hash_out = 0;

  if (b->core.n_cigar == 0 || (b->core.flag & BAM_FUNMAP)) {
    hts_pos_t pos = b->core.pos;
    hts_pos_t end_pos = bam_endpos(b);
    hts_pos_t sbin = (pos / locus_bin) * locus_bin;
    hts_pos_t ebin = (end_pos / locus_bin) * locus_bin;
    uint64_t h = 1469598103934665603ull;
    h = hash_u64pair(h, (uint64_t)sbin, (uint64_t)ebin);
    *locus_hash_out = h;
    return 0;
  }

  uint32_t* cigar = bam_get_cigar(b);
  hts_pos_t refpos = b->core.pos;

  uint64_t h = 1469598103934665603ull;
  int has_N = 0;

  for (int i=0; i < (int)b->core.n_cigar; ++i) {
    int op  = bam_cigar_op(cigar[i]);
    uint32_t len = bam_cigar_oplen(cigar[i]);
    if (op == BAM_CREF_SKIP) {
      has_N = 1;
      hts_pos_t start = refpos;
      hts_pos_t end   = refpos + len;
      if (sj_jitter>1){ start = round_by_jitter(start, sj_jitter); end = round_by_jitter(end, sj_jitter); }
      h = hash_u64pair(h, (uint64_t)start, (uint64_t)end);
      refpos += len;
    } else if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
      refpos += len;
    }
  }

  if (has_N) {
    *is_spliced = 1;
    *sj_hash_out = h;
  } else {
    hts_pos_t pos = b->core.pos;
    hts_pos_t end_pos = bam_endpos(b);
    hts_pos_t sbin = (pos / locus_bin) * locus_bin;
    hts_pos_t ebin = (end_pos / locus_bin) * locus_bin;
    uint64_t hx = 1469598103934665603ull;
    hx = hash_u64pair(hx, (uint64_t)sbin, (uint64_t)ebin);
    *locus_hash_out = hx;
  }
  return 0;
}
