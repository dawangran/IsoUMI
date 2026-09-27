#include "key_build.h"
#include "sj.h"
#include "key.h"
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static hts_pos_t bin_coord(hts_pos_t x, int bin){ return bin > 0 ? (x / bin) * bin : x; }

char* escape_key_component(const char* value){
  static const char hex[] = "0123456789ABCDEF";
  size_t need = 0;
  char* out;
  char* dst;
  if (!value) return NULL;
  for (const unsigned char* p=(const unsigned char*)value; *p; ++p){
    size_t add = (*p == '%' || *p == '|' || *p == '=') ? 3u : 1u;
    if (need > SIZE_MAX - add) return NULL;
    need += add;
  }
  /* Keep the reserved missing-value marker distinct from a literal tag value. */
  if (strcmp(value, "NA") == 0) need = 6;
  out = (char*)malloc(need + 1);
  if (!out) return NULL;
  dst = out;
  if (strcmp(value, "NA") == 0){
    memcpy(dst, "%4E%41", 6);
    dst += 6;
  } else {
    for (const unsigned char* p=(const unsigned char*)value; *p; ++p){
      if (*p == '%' || *p == '|' || *p == '='){
        *dst++ = '%';
        *dst++ = hex[*p >> 4];
        *dst++ = hex[*p & 15];
      } else {
        *dst++ = (char)*p;
      }
    }
  }
  *dst = '\0';
  return out;
}

char* build_group_key(bam1_t* b, const char* cell_tag, const char* gene_tag,
                      const char* source_tag, const char* input_scope_tag,
                      int no_gene, int no_structure, int locus_bin, int sj_jitter, int end_bin,
                      strand_mode_t strand_mode){
  const char* cb = get_tag_Z(b, cell_tag);
  const uint8_t* gene_aux = (!no_gene && gene_tag) ? bam_aux_get(b, gene_tag) : NULL;
  const char* gx = gene_aux ? bam_aux2Z(gene_aux) : NULL;
  const char* src = source_tag ? get_tag_Z(b, source_tag) : NULL;
  const char* input_scope = input_scope_tag ? get_tag_Z(b, input_scope_tag) : NULL;
  int tid = b ? b->core.tid : -1;
  int is_spliced=0; uint64_t sjh=0, lxh=0; build_sj_or_locus(b, locus_bin, sj_jitter, &is_spliced, &sjh, &lxh);
  char* cb_escaped = cb ? escape_key_component(cb) : NULL;
  char* gx_escaped = (!no_gene && gx) ? escape_key_component(gx) : NULL;
  char* src_escaped = (source_tag && src) ? escape_key_component(src) : NULL;
  char* input_escaped = (input_scope_tag && input_scope) ? escape_key_component(input_scope) : NULL;
  const char* gxv = (!no_gene && gx) ? gx_escaped : "NA";
  const char* cbv = cb ? cb_escaped : "NA";
  const char* src_prefix = source_tag ? "|SRC=" : "";
  const char* srcv = source_tag ? (src ? src_escaped : "NA") : "";
  const char* input_prefix = input_scope_tag ? "|IN=" : "";
  const char* inputv = input_scope_tag ? (input_scope ? input_escaped : "NA") : "";
  /* BAM_FREVERSE describes the alignment, not the originating RNA strand.
     A nonempty gene assignment already scopes auto mode to an annotated gene.
     Without that scope, retain the legacy separation unless explicitly disabled. */
  int has_gene_assignment = gene_aux && *gene_aux == 'Z' && gx && gx[0];
  int use_alignment_strand = strand_mode == STRAND_ALIGNMENT ||
    (strand_mode == STRAND_AUTO && !has_gene_assignment);
  int reverse = use_alignment_strand && b && (b->core.flag & BAM_FREVERSE);
  char strch = use_alignment_strand ? (reverse ? '-' : '+') : '.';
  const char* key_fmt = no_structure
    ? "CB=%s|GX=%s%s%s%s%s|TID=%d|STR=%c|CTX=NA"
    : end_bin > 0
    ? (use_alignment_strand
       ? "CB=%s|GX=%s%s%s%s%s|TID=%d|STR=%c|%s=%016llx|E5=%" PRId64 "|E3=%" PRId64
       : "CB=%s|GX=%s%s%s%s%s|TID=%d|STR=%c|%s=%016llx|EL=%" PRId64 "|ER=%" PRId64)
    : "CB=%s|GX=%s%s%s%s%s|TID=%d|STR=%c|%s=%016llx";
  hts_pos_t pos = b ? b->core.pos : -1;
  hts_pos_t end_pos = b ? bam_endpos(b) : -1;
  hts_pos_t tx5 = reverse ? end_pos : pos;
  hts_pos_t tx3 = reverse ? pos : end_pos;
  hts_pos_t e5 = bin_coord(tx5, end_bin);
  hts_pos_t e3 = bin_coord(tx3, end_bin);
  if ((cb && !cb_escaped) || (!no_gene && gx && !gx_escaped) ||
      (source_tag && src && !src_escaped) ||
      (input_scope_tag && input_scope && !input_escaped)){
    free(cb_escaped); free(gx_escaped); free(src_escaped); free(input_escaped);
    return NULL;
  }
  int need = no_structure
    ? snprintf(NULL, 0, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch)
    : end_bin > 0
    ? snprintf(NULL, 0, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch,
               is_spliced ? "SJ" : "LX",
               (unsigned long long)(is_spliced ? sjh : lxh), e5, e3)
    : snprintf(NULL, 0, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch,
               is_spliced ? "SJ" : "LX",
               (unsigned long long)(is_spliced ? sjh : lxh));
  if (need < 0){
    free(cb_escaped); free(gx_escaped); free(src_escaped); free(input_escaped);
    return NULL;
  }
  char* out = (char*)malloc((size_t)need + 1);
  if (!out){
    free(cb_escaped); free(gx_escaped); free(src_escaped); free(input_escaped);
    return NULL;
  }
  if (no_structure){
    snprintf(out, (size_t)need + 1, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch);
  } else if (end_bin > 0){
    snprintf(out, (size_t)need + 1, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch,
             is_spliced ? "SJ" : "LX",
             (unsigned long long)(is_spliced ? sjh : lxh), e5, e3);
  } else {
    snprintf(out, (size_t)need + 1, key_fmt, cbv, gxv, src_prefix, srcv, input_prefix, inputv, tid, strch,
             is_spliced ? "SJ" : "LX",
             (unsigned long long)(is_spliced ? sjh : lxh));
  }
  free(cb_escaped); free(gx_escaped); free(src_escaped); free(input_escaped);
  return out;
}
