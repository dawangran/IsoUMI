#include "pipeline.h"
#include "io.h"
#include "sj.h"
#include "umi.h"
#include "uf.h"
#include "key_build.h"
#include "key.h"
#include "vector.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static inline void ts_now(char* buf, size_t n){
  time_t t = time(NULL);
  struct tm* ptm = localtime(&t);
  if (ptm) strftime(buf, n, "%Y-%m-%d %H:%M:%S", ptm);
  else snprintf(buf, n, "0000-00-00 00:00:00");
}
#define LOG_PIPELINE(fmt, ...) do { char __ts[32]; ts_now(__ts, sizeof(__ts)); fprintf(stderr, "[%s] " fmt "\n", __ts, ##__VA_ARGS__); } while(0)
#define LOG_BUCKET(bid, fmt, ...) do { char __ts[32]; ts_now(__ts, sizeof(__ts)); fprintf(stderr, "[%s] [bucket %05d] " fmt "\n", __ts, (int)(bid), ##__VA_ARGS__); } while(0)

static int is_unmapped(const bam1_t* b){
  return !b || (b->core.flag & BAM_FUNMAP) || b->core.tid < 0;
}

static int mark_header_sort_unknown(bam_hdr_t* hdr){
  if (!hdr) return -1;
  if (sam_hdr_update_hd(hdr, "SO", "unknown") == 0) return 0;
  if (sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unknown", NULL) == 0) return 0;
  LOG_PIPELINE("failed to update output header sort order");
  return -1;
}

static void close_writers(samFile** outs, int n){
  if (!outs) return;
  for (int i=0; i<n; ++i){
    if (outs[i]) sam_close(outs[i]);
  }
}

static int ensure_dir(const char* d){
#ifdef _WIN32
  int r = mkdir(d);
#else
  int r = mkdir(d, 0775);
#endif
  if (r!=0 && errno!=EEXIST){ LOG_PIPELINE("mkdir %s failed", d); return -1; }
  return 0;
}
static unsigned long hash_cb(const char* s){
  unsigned long h=1469598103934665603ull; if (!s) return 0;
  for (const unsigned char* p=(const unsigned char*)s; *p; ++p){ h ^= *p; h *= 1099511628211ull; }
  return h;
}

static int headers_compatible(const bam_hdr_t* a, const bam_hdr_t* b){
  if (!a || !b) return 0;
  if (a->n_targets != b->n_targets) return 0;
  for (int i=0;i<a->n_targets;++i){
    if (a->target_len[i] != b->target_len[i]) return 0;
    if (strcmp(a->target_name[i], b->target_name[i]) != 0) return 0;
  }
  return 1;
}

static int header_has_exact_line(bam_hdr_t* hdr, const char* line, size_t len){
  const char* text = sam_hdr_str(hdr);
  const char* p = text;
  if (!text) return 0;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t n = e ? (size_t)(e - p) : strlen(p);
    if (n == len && strncmp(p, line, len) == 0) return 1;
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int line_type_is(const char* line, size_t len, const char* type){
  return len >= 3 && line[0] == '@' && line[1] == type[0] && line[2] == type[1];
}

static int extract_header_tag_value(const char* line, size_t len, const char* tag, char* out, size_t out_n){
  size_t tag_n = strlen(tag);
  const char* p = line;
  const char* end = line + len;
  if (out_n == 0) return -1;
  while (p < end){
    const char* next = memchr(p, '\t', (size_t)(end - p));
    size_t field_n = next ? (size_t)(next - p) : (size_t)(end - p);
    if (field_n > tag_n + 1 && strncmp(p, tag, tag_n) == 0 && p[tag_n] == ':'){
      size_t value_n = field_n - tag_n - 1;
      if (value_n + 1 > out_n) return -1;
      memcpy(out, p + tag_n + 1, value_n);
      out[value_n] = '\0';
      return 0;
    }
    if (!next) break;
    p = next + 1;
  }
  return -1;
}

static int merge_header_line(bam_hdr_t* dst, const char* line, size_t len){
  kstring_t existing = {0, 0, NULL};
  char id_buf[1024];
  char type[3];
  int find_rc;
  int ok = 0;

  if (len == 0 || line[0] != '@') return 0;
  if (line_type_is(line, len, "HD") || line_type_is(line, len, "SQ")) return 0;

  if (line_type_is(line, len, "RG") || line_type_is(line, len, "PG")){
    type[0] = line[1];
    type[1] = line[2];
    type[2] = '\0';
    if (extract_header_tag_value(line, len, "ID", id_buf, sizeof(id_buf)) != 0){
      LOG_PIPELINE("header line missing ID: %.*s", (int)len, line);
      return -1;
    }
    find_rc = sam_hdr_find_line_id(dst, type, "ID", id_buf, &existing);
    if (find_rc == 0){
      if (existing.l == len && strncmp(existing.s, line, len) == 0){
        ok = 0;
      } else {
        LOG_PIPELINE("conflicting @%c%c header line for ID=%s", type[0], type[1], id_buf);
        ok = -1;
      }
    } else if (find_rc == -1){
      ok = sam_hdr_add_lines(dst, line, len);
      if (ok != 0) LOG_PIPELINE("failed to append @%c%c header line for ID=%s", type[0], type[1], id_buf);
    } else {
      LOG_PIPELINE("failed to inspect @%c%c header line for ID=%s", type[0], type[1], id_buf);
      ok = -1;
    }
    free(existing.s);
    return ok;
  }

  if (header_has_exact_line(dst, line, len)) return 0;
  if (sam_hdr_add_lines(dst, line, len) != 0){
    LOG_PIPELINE("failed to append header line: %.*s", (int)len, line);
    return -1;
  }
  return 0;
}

static int merge_header_metadata(bam_hdr_t* dst, const bam_hdr_t* src){
  const char* text = sam_hdr_str((bam_hdr_t*)src);
  const char* p = text;
  if (!text) return -1;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (merge_header_line(dst, p, len) != 0) return -1;
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int build_merged_header(const cli_opts_t* o, bam_hdr_t** out_hdr){
  io_ctx_t io0 = {0};
  bam_hdr_t* merged = NULL;
  if (io_open(o->bam_list.data[0], NULL, &io0) != 0) return -1;
  merged = sam_hdr_dup(io0.hdr);
  if (!merged){ io_close(&io0); return -1; }
  io_close(&io0);

  for (size_t fi = 1; fi < o->bam_list.n; ++fi){
    io_ctx_t io = {0};
    if (io_open(o->bam_list.data[fi], NULL, &io) != 0){ bam_hdr_destroy(merged); return -1; }
    if (!headers_compatible(merged, io.hdr)){
      LOG_PIPELINE("[split] header mismatch in %s", o->bam_list.data[fi]);
      io_close(&io);
      bam_hdr_destroy(merged);
      return -1;
    }
    if (merge_header_metadata(merged, io.hdr) != 0){
      io_close(&io);
      bam_hdr_destroy(merged);
      return -1;
    }
    io_close(&io);
  }

  if (mark_header_sort_unknown(merged) != 0){
    bam_hdr_destroy(merged);
    return -1;
  }
  *out_hdr = merged;
  return 0;
}

static int phase_split(const cli_opts_t* o){
  bam_hdr_t* merged_hdr = NULL;
  if (ensure_dir(o->tmp_dir)!=0) return -1;
  if (build_merged_header(o, &merged_hdr) != 0) return -1;
  char path[4096];
  samFile** outs = (samFile**)calloc(o->buckets, sizeof(samFile*));
  if (!outs){ bam_hdr_destroy(merged_hdr); return -1; }
  for (int i=0;i<o->buckets;++i){
    snprintf(path, sizeof(path), "%s/bucket_%05d.bam", o->tmp_dir, i);
    outs[i] = sam_open(path, "wb");
    if (!outs[i]){ LOG_PIPELINE("[split] cannot create %s", path); bam_hdr_destroy(merged_hdr); close_writers(outs, i); free(outs); return -1; }
    if (sam_hdr_write(outs[i], merged_hdr) < 0){ LOG_PIPELINE("[split] hdr write fail %s", path); bam_hdr_destroy(merged_hdr); close_writers(outs, i + 1); free(outs); return -1; }
  }
  bam1_t* b = bam_init1();
  if (!b){ bam_hdr_destroy(merged_hdr); close_writers(outs, o->buckets); free(outs); return -1; }
  for (size_t fi=0; fi<o->bam_list.n; ++fi){
    const char* inbam = o->bam_list.data[fi];
    samFile* in = sam_open(inbam, "r");
    int read_rc = 0;
    int failed = 0;
    if (!in){ LOG_PIPELINE("[split] cannot open %s", inbam); bam_destroy1(b); bam_hdr_destroy(merged_hdr); close_writers(outs, o->buckets); free(outs); return -1; }
    bam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr){ LOG_PIPELINE("[split] cannot read header %s", inbam); sam_close(in); bam_destroy1(b); bam_hdr_destroy(merged_hdr); close_writers(outs, o->buckets); free(outs); return -1; }
    if (!headers_compatible(merged_hdr, hdr)){
      LOG_PIPELINE("[split] header mismatch in %s", inbam);
      bam_hdr_destroy(hdr); sam_close(in); bam_destroy1(b); bam_hdr_destroy(merged_hdr);
      close_writers(outs, o->buckets);
      free(outs);
      return -1;
    }
    while ((read_rc = sam_read1(in, hdr, b)) >= 0){
      const char* cb = get_tag_Z(b, o->cell_tag);
      unsigned long bid = hash_cb(cb) % (unsigned long)o->buckets;
      if (sam_write1(outs[bid], merged_hdr, b) < 0){ LOG_PIPELINE("[split] write fail for bucket %lu", bid); failed = 1; break; }
    }
    if (read_rc < -1){ LOG_PIPELINE("[split] read fail %s", inbam); failed = 1; }
    bam_hdr_destroy(hdr); sam_close(in);
    if (failed){
      bam_destroy1(b);
      bam_hdr_destroy(merged_hdr);
      close_writers(outs, o->buckets);
      free(outs);
      return -1;
    }
  }
  bam_destroy1(b); bam_hdr_destroy(merged_hdr);
  close_writers(outs, o->buckets);
  free(outs);
  return 0;
}

typedef struct { char* key; char* umi; char* umi_qual; } pair_t;
static int cmp_pair(const void* a, const void* b){ const pair_t* x=(const pair_t*)a; const pair_t* y=(const pair_t*)b; int c=strcmp(x->key, y->key); if(c) return c; return strcmp(x->umi, y->umi); }
typedef struct {
  char* umi;
  int count;
  double mean_qual;
  double* pos_qsum;
  int* pos_qn;
  int umi_len;
  int has_qual;
} umi_stat_t;
static int cmp_umi_stat_desc(const void* a, const void* b){
  const umi_stat_t* x=(const umi_stat_t*)a;
  const umi_stat_t* y=(const umi_stat_t*)b;
  if (x->count!=y->count) return (y->count-x->count);
  return strcmp(x->umi, y->umi);
}
typedef struct {
  char* raw_umi;
  char* corr_umi;
  int raw_count;
  int corr_count;
  int hamming;
  double raw_avgq;
  double corr_avgq;
  double confidence;
  int quality_supported;
} map_item_t;
static int cmp_map_item(const void* a, const void* b){ const map_item_t* x=(const map_item_t*)a; const map_item_t* y=(const map_item_t*)b; int c=strcmp(x->raw_umi, y->raw_umi); if(c) return c; return strcmp(x->corr_umi, y->corr_umi); }
typedef struct { char* key; map_item_t* items; int n; } key_map_t;
static int cmp_key_map(const void* a, const void* b){ const key_map_t* x=(const key_map_t*)a; const key_map_t* y=(const key_map_t*)b; return strcmp(x->key, y->key); }
static char* sdup(const char* s){ size_t n=s?strlen(s):0; char* p=(char*)malloc(n+1); if(!p) return NULL; memcpy(p,s?s:"",n+1); return p; }
static void free_pairs(pair_t* arr, size_t n){ for (size_t i=0;i<n;++i){ free(arr[i].key); free(arr[i].umi); free(arr[i].umi_qual); } free(arr); }
static void free_umi_stats(umi_stat_t* stats, int n){
  for (int i=0;i<n;++i){
    free(stats[i].pos_qsum);
    free(stats[i].pos_qn);
  }
  free(stats);
}
static void free_map_items(map_item_t* items, int n){ for (int i=0;i<n;++i){ free(items[i].raw_umi); free(items[i].corr_umi); } free(items); }
static void free_key_maps(key_map_t* km, int km_n){
  if (!km) return;
  for (int t=0;t<km_n;++t){
    free(km[t].key);
    free_map_items(km[t].items, km[t].n);
  }
  free(km);
}

static double clamp01(double x){
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

static int want_umi_quality(const cli_opts_t* o){
  return o->quality_aware || o->emit_explain || o->min_merge_confidence > 0.0;
}

static int add_umi_quality(umi_stat_t* st, const char* qual){
  int i;
  if (!qual) return 0;
  if (!st->umi) return 0;
  if (st->umi_len == 0) st->umi_len = (int)strlen(st->umi);
  if ((int)strlen(qual) != st->umi_len) return 0;
  if (!st->pos_qsum){
    st->pos_qsum = (double*)calloc((size_t)st->umi_len, sizeof(double));
    st->pos_qn = (int*)calloc((size_t)st->umi_len, sizeof(int));
    if (!st->pos_qsum || !st->pos_qn){
      free(st->pos_qsum);
      free(st->pos_qn);
      st->pos_qsum = NULL;
      st->pos_qn = NULL;
      return -1;
    }
  }
  for (i=0; i<st->umi_len; ++i){
    int q = (int)((unsigned char)qual[i]) - 33;
    st->pos_qsum[i] += (double)q;
    st->pos_qn[i] += 1;
  }
  st->has_qual = 1;
  return 0;
}

static void finalize_umi_quality(umi_stat_t* st){
  int i;
  double sum = 0.0;
  int n = 0;
  if (!st->has_qual || !st->pos_qsum || !st->pos_qn){
    st->mean_qual = -1.0;
    return;
  }
  for (i=0; i<st->umi_len; ++i){
    sum += st->pos_qsum[i];
    n += st->pos_qn[i];
  }
  st->mean_qual = n > 0 ? sum / (double)n : -1.0;
}

static double pos_mean_qual(const umi_stat_t* st, int pos){
  if (!st->has_qual || !st->pos_qsum || !st->pos_qn) return -1.0;
  if (pos < 0 || pos >= st->umi_len || st->pos_qn[pos] <= 0) return -1.0;
  return st->pos_qsum[pos] / (double)st->pos_qn[pos];
}

static int compute_umi_distance_metrics(const umi_stat_t* a, const umi_stat_t* b, int max_ham,
                                        int* ham_out, double* a_mismatch_q, double* b_mismatch_q,
                                        int* mismatch_q_n){
  size_t na, nb, i;
  int ham = 0, qn = 0;
  double aq = 0.0, bq = 0.0;
  if (!a || !b || !a->umi || !b->umi) return -1;
  na = strlen(a->umi);
  nb = strlen(b->umi);
  if (na != nb) return -1;
  for (i=0; i<na; ++i){
    if (a->umi[i] != b->umi[i]){
      double pa = pos_mean_qual(a, (int)i);
      double pb = pos_mean_qual(b, (int)i);
      ham++;
      if (max_ham >= 0 && ham > max_ham) return -1;
      if (pa >= 0.0 && pb >= 0.0){
        aq += pa;
        bq += pb;
        qn++;
      }
    }
  }
  *ham_out = ham;
  *a_mismatch_q = qn > 0 ? aq / (double)qn : -1.0;
  *b_mismatch_q = qn > 0 ? bq / (double)qn : -1.0;
  *mismatch_q_n = qn;
  return 0;
}

static double compute_merge_confidence(const cli_opts_t* o, const umi_stat_t* larger, const umi_stat_t* smaller,
                                       int ham, double larger_mq, double smaller_mq, int mismatch_q_n){
  double ratio = larger->count > 0 ? (double)smaller->count / (double)larger->count : 1.0;
  double size_component = clamp01(1.0 - ratio);
  double dist_component = o->ham > 0 ? clamp01(1.0 - ((double)ham / (double)(o->ham + 1))) : (ham == 0 ? 1.0 : 0.0);
  double qual_component = 0.5;
  if (o->quality_aware){
    if (mismatch_q_n > 0 && larger_mq >= 0.0 && smaller_mq >= 0.0){
      qual_component = clamp01((larger_mq - smaller_mq + 20.0) / 40.0);
    } else if (larger->mean_qual >= 0.0 && smaller->mean_qual >= 0.0){
      qual_component = clamp01((larger->mean_qual - smaller->mean_qual + 20.0) / 40.0);
    }
  }
  return clamp01(0.55 * size_component + 0.20 * dist_component + 0.25 * qual_component);
}

static int better_representative(const cli_opts_t* o, const umi_stat_t* candidate, const umi_stat_t* current){
  if (!current) return 1;
  if (candidate->count != current->count) return candidate->count > current->count;
  if (o->quality_aware){
    if (candidate->mean_qual > current->mean_qual) return 1;
    if (candidate->mean_qual < current->mean_qual) return 0;
  }
  return strcmp(candidate->umi, current->umi) < 0;
}

static const char* merge_reason(const map_item_t* item){
  if (!item) return "unknown";
  if (item->hamming == 0 || strcmp(item->raw_umi, item->corr_umi) == 0) return "self";
  if (item->quality_supported) return "count+quality";
  return "count+distance";
}

static int emit_correction_row(FILE* fp, const char* key, const map_item_t* item, int bucket){
  if (!fp || !item) return 0;
  if (fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.4f\t%s\t%d\n",
              key, item->raw_umi, item->corr_umi, item->raw_count, item->corr_count,
              item->hamming, item->raw_avgq, item->corr_avgq, item->confidence,
              merge_reason(item), bucket) < 0){
    return -1;
  }
  return 0;
}

static int build_bucket_mapping(const cli_opts_t* o, const char* bucket_bam, FILE* explain_fp, int bucket_id,
                                key_map_t** out_maps, int* out_nm){
  io_ctx_t io; if (io_open(bucket_bam, NULL, &io) != 0) return -1;
  bam1_t* b = bam_init1(); size_t cap=1<<18, n=0; pair_t* arr=(pair_t*)malloc(sizeof(pair_t)*cap); if(!b || !arr){ bam_destroy1(b); free(arr); io_close(&io); return -1; }
  int use_qual = want_umi_quality(o);
  int read_rc = 0;
  while ((read_rc = sam_read1(io.in, io.hdr, b)) >= 0){
    const char* umi_qual = NULL;
    if (is_unmapped(b)) continue;
    const char* umi = get_tag_Z(b, o->umi_tag); if(!umi) continue;
    const char* cb = get_tag_Z(b, o->cell_tag); if(!cb) continue;
    char* key = build_group_key(b, o->cell_tag, o->gene_tag, o->no_gene, o->locus_bin, o->sj_jitter, o->end_bin); if(!key) continue;
    if (n==cap){ size_t nc=cap*2; pair_t* na=(pair_t*)realloc(arr, sizeof(pair_t)*nc); if(!na){ free_pairs(arr,n); bam_destroy1(b); io_close(&io); return -1; } arr=na; cap=nc; }
    if (use_qual) umi_qual = get_tag_Z(b, o->umi_qual_tag);
    arr[n].key = key;
    arr[n].umi = sdup(umi);
    arr[n].umi_qual = (use_qual && umi_qual) ? sdup(umi_qual) : NULL;
    if (!arr[n].umi || ((use_qual && umi_qual) && !arr[n].umi_qual)){
      free(arr[n].key);
      free(arr[n].umi);
      free(arr[n].umi_qual);
      free_pairs(arr,n);
      bam_destroy1(b);
      io_close(&io);
      return -1;
    }
    n++;
  }
  if (read_rc < -1){ free_pairs(arr,n); bam_destroy1(b); io_close(&io); return -1; }
  bam_destroy1(b); io_close(&io);
  qsort(arr, n, sizeof(pair_t), cmp_pair);
  int km_cap= (n>0? (int)(n/8+16) : 16), km_n=0; key_map_t* km=(key_map_t*)malloc(sizeof(key_map_t)*km_cap); if(!km){ free_pairs(arr,n); return -1; }
  size_t i=0; while (i<n){
    size_t j=i+1; while (j<n && strcmp(arr[j].key, arr[i].key)==0) j++;
    int uc_cap=16, uc_n=0; umi_stat_t* uc=(umi_stat_t*)calloc((size_t)uc_cap, sizeof(umi_stat_t)); if(!uc){ free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
    size_t k=i; while (k<j){
      size_t t=k+1; while (t<j && strcmp(arr[t].umi, arr[k].umi)==0) t++;
      if (uc_n==uc_cap){
        int nc=uc_cap*2;
        umi_stat_t* nu=(umi_stat_t*)realloc(uc,sizeof(umi_stat_t)*(size_t)nc);
        if(!nu){ free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
        memset(nu + uc_cap, 0, sizeof(umi_stat_t) * (size_t)(nc - uc_cap));
        uc=nu; uc_cap=nc;
      }
      uc[uc_n].umi = arr[k].umi;
      uc[uc_n].count = (int)(t-k);
      uc[uc_n].umi_len = (int)strlen(arr[k].umi);
      while (k<t){
        if (add_umi_quality(&uc[uc_n], arr[k].umi_qual) != 0){ free_umi_stats(uc, uc_n + 1); free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
        k++;
      }
      finalize_umi_quality(&uc[uc_n]);
      uc_n++;
    }
    qsort(uc, uc_n, sizeof(umi_stat_t), cmp_umi_stat_desc);
    uf_t uf; if (uf_init(&uf, uc_n) != 0){ free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
    for (int a=0; a<uc_n; ++a){
      for (int b_i=a+1; b_i<uc_n; ++b_i){
        int ham = 0, mismatch_q_n = 0;
        double larger_mq = -1.0, smaller_mq = -1.0;
        double conf;
        int larger=uc[a].count, smaller=uc[b_i].count;
        if (larger<=0) continue;
        if ((double)smaller/(double)larger > o->ratio) continue;
        if (compute_umi_distance_metrics(&uc[a], &uc[b_i], o->ham, &ham, &larger_mq, &smaller_mq, &mismatch_q_n) != 0) continue;
        conf = compute_merge_confidence(o, &uc[a], &uc[b_i], ham, larger_mq, smaller_mq, mismatch_q_n);
        if (o->min_merge_confidence > 0.0 && conf < o->min_merge_confidence) continue;
        uf_union(&uf, a, b_i);
      }
    }
    int* rep=(int*)malloc(sizeof(int)*uc_n); if(!rep){ uf_free(&uf); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; } for(int idx=0; idx<uc_n; ++idx) rep[idx]=-1;
    for (int idx=0; idx<uc_n; ++idx){ int root=uf_find(&uf, idx);
      if (better_representative(o, &uc[idx], rep[root] == -1 ? NULL : &uc[rep[root]])){ rep[root]=idx; } }
    map_item_t* items=(map_item_t*)calloc((size_t)uc_n, sizeof(map_item_t)); int mi_n=0;
    if(!items){ uf_free(&uf); free(rep); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
    for (int idx=0; idx<uc_n; ++idx){
      int root=uf_find(&uf, idx);
      int ridx=rep[root];
      int ham = 0, mismatch_q_n = 0;
      double corr_mq = -1.0, raw_mq = -1.0;
      items[mi_n].raw_umi=sdup(uc[idx].umi);
      items[mi_n].corr_umi=sdup(uc[ridx].umi);
      if (!items[mi_n].raw_umi || !items[mi_n].corr_umi){
        free(items[mi_n].raw_umi);
        free(items[mi_n].corr_umi);
        free_map_items(items, mi_n);
        uf_free(&uf); free(rep); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1;
      }
      items[mi_n].raw_count = uc[idx].count;
      items[mi_n].corr_count = uc[ridx].count;
      items[mi_n].raw_avgq = uc[idx].mean_qual;
      items[mi_n].corr_avgq = uc[ridx].mean_qual;
      if (compute_umi_distance_metrics(&uc[ridx], &uc[idx], o->ham, &ham, &corr_mq, &raw_mq, &mismatch_q_n) != 0){
        ham = strcmp(uc[idx].umi, uc[ridx].umi) == 0 ? 0 : o->ham + 1;
        corr_mq = raw_mq = -1.0;
        mismatch_q_n = 0;
      }
      items[mi_n].hamming = ham;
      items[mi_n].quality_supported = o->quality_aware && mismatch_q_n > 0 && corr_mq >= 0.0 && raw_mq >= 0.0 && corr_mq > raw_mq;
      items[mi_n].confidence = strcmp(uc[idx].umi, uc[ridx].umi) == 0
        ? 1.0
        : compute_merge_confidence(o, &uc[ridx], &uc[idx], ham, corr_mq, raw_mq, mismatch_q_n);
      if (emit_correction_row(explain_fp, arr[i].key, &items[mi_n], bucket_id) != 0){
        free_map_items(items, mi_n + 1);
        uf_free(&uf); free(rep); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1;
      }
      mi_n++;
    }
    qsort(items, mi_n, sizeof(map_item_t), cmp_map_item);
    if (km_n==km_cap){ int nc=km_cap*2; key_map_t* nk=(key_map_t*)realloc(km,sizeof(key_map_t)*(size_t)nc); if(!nk){ free_map_items(items, mi_n); free(rep); uf_free(&uf); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; } km=nk; km_cap=nc; }
    km[km_n].key=sdup(arr[i].key);
    if (!km[km_n].key){ free_map_items(items, mi_n); free(rep); uf_free(&uf); free_umi_stats(uc, uc_n); free_pairs(arr,n); free_key_maps(km, km_n); return -1; }
    km[km_n].items=items; km[km_n].n=mi_n; km_n++;
    uf_free(&uf); free(rep); free_umi_stats(uc, uc_n); i=j;
  }
  free_pairs(arr,n);
  qsort(km, km_n, sizeof(key_map_t), cmp_key_map);
  *out_maps = km; *out_nm = km_n; return 0;
}

static const map_item_t* find_in_key_map(const key_map_t* km, int km_n, const char* key, const char* raw_umi){
  int lo=0, hi=km_n-1;
  while (lo<=hi){ int mid=(lo+hi)/2; int c=strcmp(km[mid].key, key);
    if (c==0){ int l=0, r=km[mid].n-1;
      while (l<=r){ int m=(l+r)/2; int d=strcmp(km[mid].items[m].raw_umi, raw_umi);
        if (d==0){ return &km[mid].items[m]; } if (d<0) l=m+1; else r=m-1; }
      return NULL; }
    if (c<0) lo=mid+1; else hi=mid-1; }
  return NULL;
}

typedef struct { char* key; char* corr; } seen_t;
static int seen_cmp(const void* a, const void* b){ const seen_t* x=(const seen_t*)a; const seen_t* y=(const seen_t*)b; int c=strcmp(x->key,y->key); if(c) return c; return strcmp(x->corr,y->corr); }
static int seen_find(seen_t* arr, int n, const char* key, const char* corr){
  int lo=0, hi=n-1; seen_t tmp={(char*)key,(char*)corr};
  while(lo<=hi){ int mid=(lo+hi)/2; int c=seen_cmp(&tmp,&arr[mid]); if(c==0) return mid; if(c<0) hi=mid-1; else lo=mid+1; }
  return -lo-1;
}
static void free_seen(seen_t* seen, int* counts, int seen_n){
  for (int s=0;s<seen_n;++s){ free(seen[s].key); free(seen[s].corr); }
  free(seen);
  free(counts);
}
static char* build_molecule_id(const char* cb, const char* key, const char* corr){
  const char* cbv = cb ? cb : "NA";
  const char* keyv = key ? key : "NA";
  const char* corrv = corr ? corr : "NA";
  int need = snprintf(NULL, 0, "%s|%s|%s", cbv, keyv, corrv);
  if (need < 0) return NULL;
  char* out = (char*)malloc((size_t)need + 1);
  if (!out) return NULL;
  snprintf(out, (size_t)need + 1, "%s|%s|%s", cbv, keyv, corrv);
  return out;
}

static int grow_seen_arrays(seen_t** seen, int** counts, int seen_n, int* seen_cap){
  int nc = *seen_cap ? (*seen_cap * 2) : 1024;
  seen_t* new_seen = (seen_t*)malloc(sizeof(seen_t) * (size_t)nc);
  int* new_counts = (int*)malloc(sizeof(int) * (size_t)nc);
  if (!new_seen || !new_counts){
    free(new_seen);
    free(new_counts);
    return -1;
  }
  if (*seen){
    memcpy(new_seen, *seen, sizeof(seen_t) * (size_t)seen_n);
    free(*seen);
  }
  if (*counts){
    memcpy(new_counts, *counts, sizeof(int) * (size_t)seen_n);
    free(*counts);
  }
  *seen = new_seen;
  *counts = new_counts;
  *seen_cap = nc;
  return 0;
}

static int phase_bucket_dedup(const cli_opts_t* o){
  int status = 0;
#pragma omp parallel for num_threads(o->threads) reduction(|:status)
  for (int i=0;i<o->buckets;i++){
    char ib[4096], ob[4096], molf[4096], asnf[4096], corrf[4096];
    snprintf(ib, sizeof(ib), "%s/bucket_%05d.bam", o->tmp_dir, i);
    snprintf(ob, sizeof(ob), "%s/bucket_%05d.dedup.bam", o->tmp_dir, i);
    snprintf(molf, sizeof(molf), "%s/bucket_%05d.molecules.tsv", o->tmp_dir, i);
    snprintf(asnf, sizeof(asnf), "%s/bucket_%05d.assignments.tsv", o->tmp_dir, i);
    snprintf(corrf, sizeof(corrf), "%s/bucket_%05d.corrections.tsv", o->tmp_dir, i);

    key_map_t* km=NULL; int km_n=0;
    int bucket_failed = 0;
    FILE* f_corr = NULL;
    if (o->emit_explain){
      f_corr = fopen(corrf, "w");
      if (!f_corr){
        LOG_BUCKET(i, "cannot create corrections TSV");
        status |= 1;
        continue;
      }
      fprintf(f_corr, "key\traw_umi\tcorr_umi\traw_count\tcorr_count\thamming\traw_avgq\tcorr_avgq\tconfidence\treason\tbucket\n");
    }
    if (build_bucket_mapping(o, ib, f_corr, i, &km, &km_n)!=0){
      if (f_corr){ fclose(f_corr); remove(corrf); }
      LOG_BUCKET(i, "mapping failed");
      status |= 1;
      continue;
    }
    if (f_corr) fclose(f_corr);

    io_ctx_t io = {0};
    if (io_open(ib, ob, &io) != 0){
      LOG_BUCKET(i, "open failed");
      if (o->emit_explain){ remove(corrf); }
      free_key_maps(km, km_n);
      status |= 1;
      continue;
    }

    FILE* f_mol = NULL; FILE* f_asn = NULL;
    if (o->emit_tsv){
      f_mol=fopen(molf,"w");
      f_asn=fopen(asnf,"w");
      if (!f_mol || !f_asn){
        LOG_BUCKET(i, "cannot create TSV outputs");
        bucket_failed = 1;
      } else {
        fprintf(f_mol,"key\tumi_corr\tcount\tbucket\n");
        fprintf(f_asn,"qname\tmolecule_id\tdup\tbucket\n");
      }
    }

    seen_t* seen=NULL; int seen_n=0, seen_cap=0;
    int* counts=NULL;
    bam1_t* b = bam_init1();
    if (!b){
      if (f_mol) fclose(f_mol);
      if (f_asn) fclose(f_asn);
      if (o->emit_explain){ remove(corrf); }
      io_close(&io);
      free_key_maps(km, km_n);
      status |= 1;
      LOG_BUCKET(i, "out of memory allocating BAM record");
      continue;
    }
    int read_rc = 0;

    while (!bucket_failed && (read_rc = sam_read1(io.in, io.hdr, b)) >= 0){
      const char* raw = get_tag_Z(b, o->umi_tag);
      const char* cb  = get_tag_Z(b, o->cell_tag);
      const map_item_t* item = NULL;
      const char* corr = NULL;
      if (is_unmapped(b)) {
        if (raw && set_tag_Z(b, o->umi_out, raw) < 0){ LOG_BUCKET(i, "failed to set passthrough UMI tag"); bucket_failed = 1; break; }
        if (set_tag_i(b, o->dup_flag, 0) < 0){ LOG_BUCKET(i, "failed to set duplicate tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      if (!raw) {
        if (set_tag_i(b, o->dup_flag, 0) < 0){ LOG_BUCKET(i, "failed to set duplicate tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      if (!cb) {
        if (set_tag_Z(b, o->umi_out, raw) < 0){ LOG_BUCKET(i, "failed to set passthrough UMI tag"); bucket_failed = 1; break; }
        if (set_tag_i(b, o->dup_flag, 0) < 0){ LOG_BUCKET(i, "failed to set duplicate tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      char* cb_copy = sdup(cb);
      if (!cb_copy){ LOG_BUCKET(i, "failed to copy cell barcode"); bucket_failed = 1; break; }
      char* key = build_group_key(b, o->cell_tag, o->gene_tag, o->no_gene, o->locus_bin, o->sj_jitter, o->end_bin);
      if (!key){ LOG_BUCKET(i, "key build failed"); free(cb_copy); bucket_failed = 1; break; }

      item = find_in_key_map(km, km_n, key, raw);
      corr = item ? item->corr_umi : raw;
      if (corr && set_tag_Z(b, o->umi_out, corr) < 0){ LOG_BUCKET(i, "failed to set corrected UMI tag"); free(cb_copy); free(key); bucket_failed = 1; break; }

      int idx = seen_find(seen, seen_n, key, corr);
      int is_dup = 0;
      if (idx>=0){
        counts[idx]++;
        is_dup = 1;
      } else {
        int ins = -idx-1;
        if (seen_n==seen_cap){
          if (grow_seen_arrays(&seen, &counts, seen_n, &seen_cap) != 0){
            free(cb_copy);
            free(key);
            bucket_failed = 1;
            break;
          }
        }
        memmove(seen + ins + 1, seen + ins, sizeof(seen_t)*(seen_n - ins));
        memmove(counts + ins + 1, counts + ins, sizeof(int)*(seen_n - ins));
        seen[ins].key = sdup(key);
        seen[ins].corr = sdup(corr);
        if (!seen[ins].key || !seen[ins].corr){
          free(seen[ins].key);
          free(seen[ins].corr);
          free(cb_copy);
          free(key);
          bucket_failed = 1;
          break;
        }
        counts[ins]=1; seen_n++;
      }

      if (set_tag_i(b, o->dup_flag, is_dup ? 1 : 0) < 0){ LOG_BUCKET(i, "failed to set duplicate tag"); free(cb_copy); free(key); bucket_failed = 1; break; }
      if (o->mol_tag){
        char* mi = build_molecule_id(cb_copy, key, corr);
        if (!mi || set_tag_Z(b, o->mol_tag, mi) < 0){
          free(mi);
          free(cb_copy);
          free(key);
          LOG_BUCKET(i, "failed to build molecule id");
          bucket_failed = 1;
          break;
        }
        free(mi);
      }
      if (o->emit_tsv && f_asn){
        const char* qn = bam_get_qname(b);
        char* mi2 = build_molecule_id(cb_copy, key, corr);
        if (!mi2){
          free(cb_copy);
          free(key);
          LOG_BUCKET(i, "failed to build assignment id");
          bucket_failed = 1;
          break;
        }
        fprintf(f_asn, "%s\t%s\t%d\t%d\n", qn?qn:"*", mi2, is_dup ? 1 : 0, i);
        free(mi2);
      }

      if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); free(cb_copy); free(key); bucket_failed = 1; break; }
      free(cb_copy);
      free(key);
    }
    if (read_rc < -1){ LOG_BUCKET(i, "read fail"); bucket_failed = 1; }
    bam_destroy1(b); io_close(&io);

    if (!bucket_failed && o->emit_tsv && f_mol){
      for (int t=0;t<seen_n;++t){
        fprintf(f_mol, "%s\t%s\t%d\t%d\n", seen[t].key, seen[t].corr, counts[t], i);
      }
    }
    if (f_mol) fclose(f_mol);
    if (f_asn) fclose(f_asn);

    free_key_maps(km, km_n);
    free_seen(seen, counts, seen_n);

    if (bucket_failed){
      remove(ob);
      if (o->emit_tsv){ remove(molf); remove(asnf); }
      if (o->emit_explain){ remove(corrf); }
      status |= 1;
      LOG_BUCKET(i, "failed");
    } else {
      LOG_BUCKET(i, "done");
    }
  }
  return status ? -1 : 0;
}

static int phase_concat(const cli_opts_t* o){
  char outp[4096]; snprintf(outp, sizeof(outp), "%s.dedup.bam", o->out_prefix);
  samFile* out = sam_open(outp, "wb"); if (!out){ LOG_PIPELINE("[concat] cannot create %s", outp); return -1; }
  char hb[4096]; snprintf(hb, sizeof(hb), "%s/bucket_%05d.dedup.bam", o->tmp_dir, 0);
  samFile* hin = sam_open(hb, "r"); if (!hin){ LOG_PIPELINE("[concat] cannot open header %s", hb); sam_close(out); return -1; }
  bam_hdr_t* hdr = sam_hdr_read(hin); if (!hdr){ LOG_PIPELINE("[concat] cannot read header"); sam_close(hin); sam_close(out); return -1; }
  sam_close(hin);
  if (mark_header_sort_unknown(hdr) != 0){ bam_hdr_destroy(hdr); sam_close(out); return -1; }
  if (sam_hdr_write(out, hdr) < 0){ LOG_PIPELINE("[concat] hdr write fail"); bam_hdr_destroy(hdr); sam_close(out); return -1; }
  bam1_t* b = bam_init1();
  if (!b){ bam_hdr_destroy(hdr); sam_close(out); return -1; }
  int failed = 0;
  for (int i=0;i<o->buckets;i++){
    char ib[4096]; snprintf(ib, sizeof(ib), "%s/bucket_%05d.dedup.bam", o->tmp_dir, i);
    samFile* in = sam_open(ib, "r"); if (!in){ LOG_PIPELINE("[concat] cannot open %s", ib); failed = 1; break; }
    bam_hdr_t* h2 = sam_hdr_read(in);
    int read_rc = 0;
    if (!h2){ LOG_PIPELINE("[concat] cannot read header %s", ib); sam_close(in); failed = 1; break; }
    bam_hdr_destroy(h2);
    while ((read_rc = sam_read1(in, hdr, b)) >= 0){
      if (sam_write1(out, hdr, b) < 0){ LOG_PIPELINE("[concat] write fail"); failed = 1; break; }
    }
    if (read_rc < -1){ LOG_PIPELINE("[concat] read fail %s", ib); failed = 1; }
    sam_close(in);
    if (failed) break;
  }
  bam_destroy1(b); bam_hdr_destroy(hdr); sam_close(out); return failed ? -1 : 0;
}

static void cleanup_tmp(const cli_opts_t* o){
  for (int i=0;i<o->buckets;i++){
    char p1[4096], p2[4096], p3[4096], p4[4096], p5[4096];
    snprintf(p1, sizeof(p1), "%s/bucket_%05d.bam",         o->tmp_dir, i);
    snprintf(p2, sizeof(p2), "%s/bucket_%05d.dedup.bam",   o->tmp_dir, i);
    snprintf(p3, sizeof(p3), "%s/bucket_%05d.molecules.tsv",   o->tmp_dir, i);
    snprintf(p4, sizeof(p4), "%s/bucket_%05d.assignments.tsv", o->tmp_dir, i);
    snprintf(p5, sizeof(p5), "%s/bucket_%05d.corrections.tsv", o->tmp_dir, i);
    remove(p1); remove(p2); remove(p3); remove(p4); remove(p5);
  }
  rmdir(o->tmp_dir);
}

static int concat_bucket_report(const cli_opts_t* o, const char* suffix, const char* header_prefix, const char* output_name){
  FILE* out = NULL;
  char outp[4096];
  snprintf(outp, sizeof(outp), "%s.%s", o->out_prefix, output_name);
  out = fopen(outp, "w");
  if (!out){ LOG_PIPELINE("[concat] cannot create %s", outp); return -1; }
  for (int i=0;i<o->buckets;i++){
    char path[4096];
    FILE* in = NULL;
    char buf[1<<15];
    snprintf(path, sizeof(path), "%s/bucket_%05d.%s", o->tmp_dir, i, suffix);
    in = fopen(path, "r");
    if (!in){ LOG_PIPELINE("[concat] cannot open %s", path); fclose(out); return -1; } else {
        if (fgets(buf,sizeof(buf),in)){
          if (i == 0 || strncmp(buf, header_prefix, strlen(header_prefix)) != 0) fputs(buf,out);
        }
        while (fgets(buf,sizeof(buf),in)) fputs(buf,out);
        if (ferror(in)){ LOG_PIPELINE("[concat] read fail %s", path); fclose(in); fclose(out); return -1; }
        fclose(in);
    }
  }
  fclose(out);
  return 0;
}

static int concat_reports(const cli_opts_t* o){
  if (o->emit_tsv){
    if (concat_bucket_report(o, "molecules.tsv", "key\t", "molecules.tsv") != 0) return -1;
    if (concat_bucket_report(o, "assignments.tsv", "qname\t", "assignments.tsv") != 0) return -1;
  }
  if (o->emit_explain){
    if (concat_bucket_report(o, "corrections.tsv", "key\t", "corrections.tsv") != 0) return -1;
  }
  return 0;
}

int run_pipeline(const cli_opts_t* o){
  LOG_PIPELINE("Phase 1: split into %d buckets at %s", o->buckets, o->tmp_dir);
  if (phase_split(o)!=0) return -1;
  LOG_PIPELINE("Phase 2: per-bucket dedup (threads=%d)", o->threads);
  if (phase_bucket_dedup(o)!=0) return -1;
  LOG_PIPELINE("Phase 3: concat into %s.dedup.bam", o->out_prefix);
  if (phase_concat(o)!=0) return -1;
  if (o->emit_tsv || o->emit_explain){
    LOG_PIPELINE("Phase 4: concat reports for %s", o->out_prefix);
    if (concat_reports(o)!=0) return -1;
  }
  if (!o->keep_tmp){
    LOG_PIPELINE("Cleanup: removing temp files under %s", o->tmp_dir);
    cleanup_tmp(o);
  }
  LOG_PIPELINE("DONE");
  return 0;
}
