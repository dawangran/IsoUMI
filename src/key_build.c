#include "key_build.h"
#include "sj.h"
#include "key.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static char strand_char(bam1_t* b){ return (b->core.flag & BAM_FREVERSE) ? '-' : '+'; }
static int bin_coord(int x, int bin){ return bin > 0 ? (x / bin) * bin : x; }

char* build_group_key(bam1_t* b, const char* cell_tag, const char* gene_tag,
                      int no_gene, int locus_bin, int sj_jitter, int end_bin){
  const char* cb = get_tag_Z(b, cell_tag);
  const char* gx = NULL; if (!no_gene) gx = get_tag_Z(b, gene_tag);
  int tid = b ? b->core.tid : -1;
  int is_spliced=0; uint64_t sjh=0, lxh=0; build_sj_or_locus(b, locus_bin, sj_jitter, &is_spliced, &sjh, &lxh);
  const char* gxv = (!no_gene && gx) ? gx : "NA";
  const char* cbv = cb?cb:"NA";
  char strch = strand_char(b);
  const char* key_fmt = end_bin > 0
    ? "CB=%s|GX=%s|TID=%d|STR=%c|%s=%016llx|E5=%d|E3=%d"
    : "CB=%s|GX=%s|TID=%d|STR=%c|%s=%016llx";
  int32_t pos = b ? b->core.pos : -1;
  int32_t end_pos = b ? bam_endpos(b) : -1;
  int tx5 = (b && (b->core.flag & BAM_FREVERSE)) ? end_pos : pos;
  int tx3 = (b && (b->core.flag & BAM_FREVERSE)) ? pos : end_pos;
  int e5 = bin_coord(tx5, end_bin);
  int e3 = bin_coord(tx3, end_bin);
  int need = end_bin > 0
    ? snprintf(NULL, 0, key_fmt, cbv, gxv, tid, strch,
               is_spliced ? "SJ" : "LX",
               (unsigned long long)(is_spliced ? sjh : lxh), e5, e3)
    : snprintf(NULL, 0, key_fmt, cbv, gxv, tid, strch,
               is_spliced ? "SJ" : "LX",
               (unsigned long long)(is_spliced ? sjh : lxh));
  if (need < 0) return NULL;
  char* out = (char*)malloc((size_t)need + 1);
  if (!out) return NULL;
  if (end_bin > 0){
    snprintf(out, (size_t)need + 1, key_fmt, cbv, gxv, tid, strch,
             is_spliced ? "SJ" : "LX",
             (unsigned long long)(is_spliced ? sjh : lxh), e5, e3);
  } else {
    snprintf(out, (size_t)need + 1, key_fmt, cbv, gxv, tid, strch,
             is_spliced ? "SJ" : "LX",
             (unsigned long long)(is_spliced ? sjh : lxh));
  }
  return out;
}
