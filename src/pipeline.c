#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200809L
#endif

#include "pipeline.h"
#include "io.h"
#include "sj.h"
#include "umi.h"
#include "key_build.h"
#include "key.h"
#include "vector.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <htslib/kstring.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static inline void ts_now(char* buf, size_t n){
  time_t t = time(NULL);
  struct tm tm_value;
#ifdef _WIN32
  struct tm* ptm = localtime_s(&tm_value, &t) == 0 ? &tm_value : NULL;
#else
  struct tm* ptm = localtime_r(&t, &tm_value);
#endif
  if (ptm) strftime(buf, n, "%Y-%m-%d %H:%M:%S", ptm);
  else snprintf(buf, n, "0000-00-00 00:00:00");
}
#define LOG_PIPELINE(fmt, ...) do { char __ts[32]; ts_now(__ts, sizeof(__ts)); fprintf(stderr, "[%s] " fmt "\n", __ts, ##__VA_ARGS__); } while(0)
#define LOG_BUCKET(bid, fmt, ...) do { char __ts[32]; ts_now(__ts, sizeof(__ts)); fprintf(stderr, "[%s] [bucket %05d] " fmt "\n", __ts, (int)(bid), ##__VA_ARGS__); } while(0)

static int is_unmapped(const bam1_t* b){
  return !b || (b->core.flag & BAM_FUNMAP) || b->core.tid < 0;
}

static int is_non_primary(const bam1_t* b){
  return b && (b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY));
}

static int mark_header_sort_unknown(bam_hdr_t* hdr){
  if (!hdr) return -1;
  if (sam_hdr_update_hd(hdr, "SO", "unknown") == 0) return 0;
  if (sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unknown", NULL) == 0) return 0;
  LOG_PIPELINE("failed to update output header sort order");
  return -1;
}

static int close_writers(samFile** outs, int n){
  int failed = 0;
  if (!outs) return 0;
  for (int i=0; i<n; ++i){
    if (outs[i] && sam_close(outs[i]) < 0) failed = 1;
    outs[i] = NULL;
  }
  return failed ? -1 : 0;
}

static int ensure_dir(const char* d){
#ifdef _WIN32
  int r = mkdir(d);
#else
  int r = mkdir(d, 0775);
#endif
  if (r!=0){
    if (errno==EEXIST) LOG_PIPELINE("temporary directory already exists: %s", d);
    else LOG_PIPELINE("mkdir %s failed", d);
    return -1;
  }
  return 0;
}

static char* path_with_suffix(const char* prefix, const char* suffix){
  size_t a;
  size_t b;
  char* out;
  if (!prefix || !suffix) return NULL;
  a = strlen(prefix);
  b = strlen(suffix);
  if (a > SIZE_MAX - b - 1) return NULL;
  out = (char*)malloc(a + b + 1);
  if (!out) return NULL;
  memcpy(out, prefix, a);
  memcpy(out + a, suffix, b + 1);
  return out;
}

static int paths_refer_to_same_file(const char* a, const char* b){
  struct stat sa;
  struct stat sb;
  if (!a || !b) return 0;
  if (strcmp(a, b) == 0) return 1;
  if (stat(a, &sa) != 0 || stat(b, &sb) != 0) return 0;
  return sa.st_dev == sb.st_dev && sa.st_ino == sb.st_ino;
}

static int validate_output_aliases(const cli_opts_t* o){
  const char* suffixes[4];
  size_t suffix_n = 0;
  suffixes[suffix_n++] = ".dedup.bam";
  if (o->emit_tsv){
    suffixes[suffix_n++] = ".molecules.tsv";
    suffixes[suffix_n++] = ".assignments.tsv";
  }
  if (o->emit_explain) suffixes[suffix_n++] = ".corrections.tsv";
  for (size_t si=0; si<suffix_n; ++si){
    char* output_path = path_with_suffix(o->out_prefix, suffixes[si]);
    if (!output_path){
      LOG_PIPELINE("cannot allocate output path");
      return -1;
    }
    for (size_t fi=0; fi<o->bam_list.n; ++fi){
      if (paths_refer_to_same_file(output_path, o->bam_list.data[fi])){
        LOG_PIPELINE("refusing to overwrite input %s with output %s",
                     o->bam_list.data[fi], output_path);
        free(output_path);
        return -1;
      }
    }
    free(output_path);
  }
  return 0;
}

static int create_staging_path(const char* final_path, char** out_path){
  size_t base_len;
  char* path;
  int fd;
  if (!final_path || !out_path) return -1;
  base_len = strlen(final_path);
  if (base_len > SIZE_MAX - 64) return -1;
  path = (char*)malloc(base_len + 64);
  if (!path) return -1;
  for (unsigned int attempt=0; attempt<1000; ++attempt){
    int n = snprintf(path, base_len + 64, "%s.partial.%ld.%u",
                     final_path, (long)getpid(), attempt);
    if (n < 0 || (size_t)n >= base_len + 64){
      free(path);
      return -1;
    }
    fd = open(path, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd >= 0){
      if (close(fd) != 0){
        remove(path);
        free(path);
        return -1;
      }
      *out_path = path;
      return 0;
    }
    if (errno != EEXIST){
      LOG_PIPELINE("cannot create staging file %s: %s", path, strerror(errno));
      free(path);
      return -1;
    }
  }
  LOG_PIPELINE("cannot allocate a unique staging path for %s", final_path);
  free(path);
  return -1;
}

static unsigned long hash_cb(const char* s){
  unsigned long h=1469598103934665603ull; if (!s) return 0;
  for (const unsigned char* p=(const unsigned char*)s; *p; ++p){ h ^= *p; h *= 1099511628211ull; }
  return h;
}

typedef struct {
  char* old_id;
  char* new_id;
} pg_id_map_item_t;

typedef struct {
  pg_id_map_item_t* items;
  size_t n;
  size_t cap;
} pg_id_map_t;

static char* dup_slice(const char* s, size_t n){
  char* out = (char*)malloc(n + 1);
  if (!out) return NULL;
  memcpy(out, s, n);
  out[n] = '\0';
  return out;
}
typedef struct {
  char** slots;
  size_t n;
  size_t cap;
} string_set_t;

static void string_set_destroy(string_set_t* set){
  if (!set) return;
  for (size_t i=0; i<set->cap; ++i) free(set->slots[i]);
  free(set->slots);
  set->slots = NULL;
  set->n = set->cap = 0;
}

static int string_set_rehash(string_set_t* set, size_t new_cap){
  char** slots;
  if (!set || new_cap < 16 || (new_cap & (new_cap - 1)) != 0) return -1;
  if (new_cap > SIZE_MAX / sizeof(*slots)) return -1;
  slots = (char**)calloc(new_cap, sizeof(*slots));
  if (!slots) return -1;
  for (size_t i=0; i<set->cap; ++i){
    char* value = set->slots[i];
    if (value){
      size_t idx = (size_t)hash_cb(value) & (new_cap - 1);
      while (slots[idx]) idx = (idx + 1) & (new_cap - 1);
      slots[idx] = value;
    }
  }
  free(set->slots);
  set->slots = slots;
  set->cap = new_cap;
  return 0;
}

static int string_set_add(string_set_t* set, const char* value){
  size_t idx;
  char* copy;
  if (!set || !value) return -1;
  if (set->cap == 0 && string_set_rehash(set, 1024) != 0) return -1;
  if (set->n >= set->cap - set->cap / 3){
    if (set->cap > SIZE_MAX / 2 ||
        string_set_rehash(set, set->cap * 2) != 0) return -1;
  }
  idx = (size_t)hash_cb(value) & (set->cap - 1);
  while (set->slots[idx]){
    if (strcmp(set->slots[idx], value) == 0) return 0;
    idx = (idx + 1) & (set->cap - 1);
  }
  copy = dup_slice(value, strlen(value));
  if (!copy) return -1;
  set->slots[idx] = copy;
  set->n++;
  return 0;
}

static int string_set_export(string_set_t* set, strvec_t* out){
  size_t j = 0;
  char** data;
  if (!set || !out || out->n != 0 || out->data) return -1;
  if (set->n == 0) return 0;
  if (set->n > SIZE_MAX / sizeof(*data)) return -1;
  data = (char**)malloc(set->n * sizeof(*data));
  if (!data) return -1;
  for (size_t i=0; i<set->cap; ++i){
    if (set->slots[i]){
      data[j++] = set->slots[i];
      set->slots[i] = NULL;
    }
  }
  out->data = data;
  out->n = out->cap = j;
  return 0;
}


static void free_pg_id_map(pg_id_map_t* map){
  if (!map) return;
  for (size_t i=0; i<map->n; ++i){
    free(map->items[i].old_id);
    free(map->items[i].new_id);
  }
  free(map->items);
  map->items = NULL;
  map->n = map->cap = 0;
}

static void free_pg_id_maps(pg_id_map_t* maps, size_t n){
  if (!maps) return;
  for (size_t i=0; i<n; ++i) free_pg_id_map(&maps[i]);
  free(maps);
}

static const char* pg_id_map_lookup(const pg_id_map_t* map, const char* old_id){
  if (!map || !old_id) return NULL;
  for (size_t i=0; i<map->n; ++i){
    if (strcmp(map->items[i].old_id, old_id) == 0) return map->items[i].new_id;
  }
  return NULL;
}

static int pg_id_map_add(pg_id_map_t* map, const char* old_id, const char* new_id){
  pg_id_map_item_t* grown;
  if (!map || !old_id || !new_id) return -1;
  if (pg_id_map_lookup(map, old_id)) return -1;
  if (map->n == map->cap){
    size_t new_cap = map->cap ? map->cap * 2 : 8;
    grown = (pg_id_map_item_t*)realloc(map->items, sizeof(*grown) * new_cap);
    if (!grown) return -1;
    map->items = grown;
    map->cap = new_cap;
  }
  map->items[map->n].old_id = dup_slice(old_id, strlen(old_id));
  map->items[map->n].new_id = dup_slice(new_id, strlen(new_id));
  if (!map->items[map->n].old_id || !map->items[map->n].new_id){
    free(map->items[map->n].old_id);
    free(map->items[map->n].new_id);
    return -1;
  }
  map->n++;
  return 0;
}

static int strvec_contains_value(const strvec_t* values, const char* value){
  if (!values || !value) return 0;
  for (size_t i=0; i<values->n; ++i){
    if (strcmp(values->data[i], value) == 0) return 1;
  }
  return 0;
}

static int input_scope_label(size_t input_idx, char* buf, size_t n){
  if (!buf || n == 0) return -1;
  if (snprintf(buf, n, "input%06zu", input_idx) < 0) return -1;
  return 0;
}

static const char* active_input_scope_tag(const cli_opts_t* o){
  return (o && o->isolate_inputs) ? o->input_scope_tag : NULL;
}

static int set_input_scope_tag(const cli_opts_t* o, bam1_t* b, const char* label){
  if (!o || !o->isolate_inputs) return 0;
  if (bam_aux_get(b, o->input_scope_tag)){
    LOG_PIPELINE("input-scope tag %s already exists; choose a different --input-scope-tag", o->input_scope_tag);
    return -1;
  }
  return set_tag_Z(b, o->input_scope_tag, label);
}

static int clear_input_scope_tag(const cli_opts_t* o, bam1_t* b){
  uint8_t* aux;
  if (!o || !o->isolate_inputs) return 0;
  aux = bam_aux_get(b, o->input_scope_tag);
  if (!aux) return 0;
  return bam_aux_del(b, aux);
}

static int set_duplicate_status(const cli_opts_t* o, bam1_t* b, int is_dup){
  if (o->set_bam_dup_flag){
    if (is_dup) b->core.flag |= BAM_FDUP;
    else b->core.flag &= (uint16_t)~(uint16_t)BAM_FDUP;
  }
  return set_tag_i(b, o->dup_flag, is_dup ? 1 : 0);
}

static char* extract_header_tag_value_dup(const char* line, size_t len, const char* tag);
static int clear_tag_if_present(bam1_t* b, const char* tag){
  uint8_t* aux;
  if (!b || !tag) return 0;
  for (;;){
    errno = 0;
    aux = bam_aux_get(b, tag);
    if (!aux) return errno == EINVAL ? -1 : 0;
    if (bam_aux_del(b, aux) != 0) return -1;
  }
}

static int clear_molecule_tag(const cli_opts_t* o, bam1_t* b){
  return o->mol_tag ? clear_tag_if_present(b, o->mol_tag) : 0;
}


static char* sq_tag_value_dup(const bam_hdr_t* hdr, const char* sequence_name,
                              const char* tag){
  kstring_t line = {0, 0, NULL};
  char* value = NULL;
  if (sam_hdr_find_line_id((bam_hdr_t*)hdr, "SQ", "SN", sequence_name, &line) == 0){
    value = extract_header_tag_value_dup(line.s, line.l, tag);
  }
  free(line.s);
  return value;
}

static int matching_optional_sq_tag(const bam_hdr_t* a, const bam_hdr_t* b,
                                    const char* sequence_name, const char* tag){
  char* av = sq_tag_value_dup(a, sequence_name, tag);
  char* bv = sq_tag_value_dup(b, sequence_name, tag);
  int matches = !av || !bv || strcmp(av, bv) == 0;
  if (!matches){
    LOG_PIPELINE("reference metadata mismatch for %s: @SQ %s:%s versus %s:%s",
                 sequence_name, tag, av, tag, bv);
  }
  free(av);
  free(bv);
  return matches;
}

static int adopt_optional_sq_tags(bam_hdr_t* dst, const bam_hdr_t* src){
  static const char* tags[] = {"M5", "AS"};
  if (!dst || !src) return -1;
  for (int i=0; i<dst->n_targets; ++i){
    const char* sequence_name = src->target_name[i];
    for (size_t j=0; j<sizeof(tags)/sizeof(tags[0]); ++j){
      char* dst_value = sq_tag_value_dup(dst, sequence_name, tags[j]);
      char* src_value = sq_tag_value_dup(src, sequence_name, tags[j]);
      if (!dst_value && src_value &&
          sam_hdr_update_line(dst, "SQ", "SN", sequence_name,
                              tags[j], src_value, NULL) != 0){
        LOG_PIPELINE("failed to preserve @SQ %s:%s for %s",
                     tags[j], src_value, sequence_name);
        free(dst_value);
        free(src_value);
        return -1;
      }
      free(dst_value);
      free(src_value);
    }
  }
  return 0;
}

static int headers_compatible(const bam_hdr_t* a, const bam_hdr_t* b){
  if (!a || !b) return 0;
  if (a->n_targets != b->n_targets) return 0;
  for (int i=0;i<a->n_targets;++i){
    if (a->target_len[i] != b->target_len[i]) return 0;
    if (strcmp(a->target_name[i], b->target_name[i]) != 0) return 0;
    if (!matching_optional_sq_tag(a, b, a->target_name[i], "M5")) return 0;
    if (!matching_optional_sq_tag(a, b, a->target_name[i], "AS")) return 0;
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

static char* extract_header_tag_value_dup(const char* line, size_t len, const char* tag){
  size_t tag_n = strlen(tag);
  const char* p = line;
  const char* end = line + len;
  while (p < end){
    const char* next = memchr(p, '\t', (size_t)(end - p));
    size_t field_n = next ? (size_t)(next - p) : (size_t)(end - p);
    if (field_n > tag_n + 1 && strncmp(p, tag, tag_n) == 0 && p[tag_n] == ':'){
      size_t value_n = field_n - tag_n - 1;
      return dup_slice(p + tag_n + 1, value_n);
    }
    if (!next) break;
    p = next + 1;
  }
  return NULL;
}

static char* make_scoped_pg_id(const char* old_id, size_t input_idx,
                               const strvec_t* used_ids){
  int attempt = 0;
  for (;;){
    int need = attempt == 0
      ? snprintf(NULL, 0, "%s.isoumi.input%06zu", old_id, input_idx)
      : snprintf(NULL, 0, "%s.isoumi.input%06zu.%d", old_id, input_idx, attempt);
    char* candidate;
    if (need < 0) return NULL;
    candidate = (char*)malloc((size_t)need + 1);
    if (!candidate) return NULL;
    if (attempt == 0){
      snprintf(candidate, (size_t)need + 1, "%s.isoumi.input%06zu", old_id, input_idx);
    } else {
      snprintf(candidate, (size_t)need + 1, "%s.isoumi.input%06zu.%d", old_id, input_idx, attempt);
    }
    if (!strvec_contains_value(used_ids, candidate)) return candidate;
    free(candidate);
    attempt++;
  }
}

static char* make_unique_pg_id(const char* base, const strvec_t* used_ids){
  int attempt = 0;
  for (;;){
    int need = attempt == 0
      ? snprintf(NULL, 0, "%s", base)
      : snprintf(NULL, 0, "%s.%d", base, attempt);
    char* candidate;
    if (need < 0) return NULL;
    candidate = (char*)malloc((size_t)need + 1);
    if (!candidate) return NULL;
    if (attempt == 0) snprintf(candidate, (size_t)need + 1, "%s", base);
    else snprintf(candidate, (size_t)need + 1, "%s.%d", base, attempt);
    if (!strvec_contains_value(used_ids, candidate)) return candidate;
    free(candidate);
    attempt++;
  }
}

static int add_isoumi_program_record(bam_hdr_t* hdr, const cli_opts_t* o,
                                     strvec_t* used_ids){
  char* id = make_unique_pg_id("IsoUMI", used_ids);
  int rc;
  if (!id) return -1;
  rc = sam_hdr_add_line(hdr, "PG",
                        "ID", id,
                        "PN", "IsoUMI",
                        "VN", ISOUMI_VERSION_NUMBER,
                        "CL", o->command_line ? o->command_line : "isoumi",
                        NULL);
  if (rc != 0){
    LOG_PIPELINE("failed to append IsoUMI @PG header line");
    free(id);
    return -1;
  }
  if (strvec_push(used_ids, id) != 0){
    free(id);
    return -1;
  }
  free(id);
  return 0;
}

static int build_pg_id_map(const bam_hdr_t* hdr, size_t input_idx,
                           strvec_t* used_ids, pg_id_map_t* map){
  const char* text = sam_hdr_str((bam_hdr_t*)hdr);
  const char* p = text;
  if (!text) return -1;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (line_type_is(p, len, "PG")){
      char* old_id = extract_header_tag_value_dup(p, len, "ID");
      char* new_id = NULL;
      if (!old_id){
        LOG_PIPELINE("@PG header line missing ID in input %zu: %.*s", input_idx, (int)len, p);
        return -1;
      }
      if (pg_id_map_lookup(map, old_id)){
        LOG_PIPELINE("duplicate @PG ID=%s within input %zu", old_id, input_idx);
        free(old_id);
        return -1;
      }
      new_id = input_idx == 0
        ? dup_slice(old_id, strlen(old_id))
        : make_scoped_pg_id(old_id, input_idx, used_ids);
      if (!new_id || pg_id_map_add(map, old_id, new_id) != 0 ||
          strvec_push(used_ids, new_id) != 0){
        free(old_id);
        free(new_id);
        return -1;
      }
      if (strcmp(old_id, new_id) != 0){
        LOG_PIPELINE("renaming @PG ID=%s from input %zu to %s", old_id, input_idx, new_id);
      }
      free(old_id);
      free(new_id);
    }
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int validate_pg_parents(const bam_hdr_t* hdr, size_t input_idx,
                               const pg_id_map_t* map){
  const char* text = sam_hdr_str((bam_hdr_t*)hdr);
  const char* p = text;
  if (!text) return -1;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (line_type_is(p, len, "PG")){
      char* parent = extract_header_tag_value_dup(p, len, "PP");
      if (parent && !pg_id_map_lookup(map, parent)){
        char* id = extract_header_tag_value_dup(p, len, "ID");
        LOG_PIPELINE("@PG ID=%s in input %zu references missing parent PP=%s",
                     id ? id : "(missing)", input_idx, parent);
        free(id);
        free(parent);
        return -1;
      }
      free(parent);
    }
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int rewrite_pg_header_line(const char* line, size_t len,
                                  const pg_id_map_t* map, kstring_t* out){
  const char* p = line;
  const char* end = line + len;
  out->l = 0;
  while (p < end){
    const char* next = memchr(p, '\t', (size_t)(end - p));
    size_t field_n = next ? (size_t)(next - p) : (size_t)(end - p);
    const char* replacement = NULL;
    if (field_n > 3 && ((p[0] == 'I' && p[1] == 'D') ||
                        (p[0] == 'P' && p[1] == 'P')) && p[2] == ':'){
      char* old_value = dup_slice(p + 3, field_n - 3);
      if (!old_value) return -1;
      replacement = pg_id_map_lookup(map, old_value);
      if (!replacement){
        LOG_PIPELINE("@PG field %c%c references missing program ID=%s",
                     p[0], p[1], old_value);
        free(old_value);
        return -1;
      }
      if (kputsn(p, 3, out) < 0 || kputs(replacement, out) < 0){
        free(old_value);
        return -1;
      }
      free(old_value);
    } else if (kputsn(p, (int)field_n, out) < 0){
      return -1;
    }
    if (!next) break;
    if (kputc('\t', out) < 0) return -1;
    p = next + 1;
  }
  return 0;
}

static int rewrite_record_pg(bam1_t* b, const pg_id_map_t* map, const char* input_name){
  uint8_t* aux = bam_aux_get(b, "PG");
  const char* old_id;
  const char* new_id;
  if (!aux) return 0;
  old_id = bam_aux2Z(aux);
  if (!old_id){
    LOG_PIPELINE("record %s has a non-string PG tag in %s", bam_get_qname(b), input_name);
    return -1;
  }
  new_id = pg_id_map_lookup(map, old_id);
  if (!new_id){
    LOG_PIPELINE("record %s references missing @PG ID=%s in %s", bam_get_qname(b), old_id, input_name);
    return -1;
  }
  if (strcmp(old_id, new_id) == 0) return 0;
  return set_tag_Z(b, "PG", new_id);
}

static int strip_record_pg(bam1_t* b, const char* input_name){
  uint8_t* aux;
  for (;;){
    errno = 0;
    aux = bam_aux_get(b, "PG");
    if (!aux){
      if (errno == EINVAL){
        LOG_PIPELINE("record %s has corrupt auxiliary data in %s",
                     bam_get_qname(b), input_name);
        return -1;
      }
      return 0;
    }
    if (bam_aux_del(b, aux) != 0){
      LOG_PIPELINE("failed to remove PG tag from record %s in %s",
                   bam_get_qname(b), input_name);
      return -1;
    }
  }
}

static bam_hdr_t* duplicate_header_without_pg(const bam_hdr_t* src){
  const char* text = sam_hdr_str((bam_hdr_t*)src);
  const char* p = text;
  kstring_t filtered = {0, 0, NULL};
  bam_hdr_t* out = NULL;
  if (!text) return NULL;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (!line_type_is(p, len, "PG")){
      if (len > INT_MAX || kputsn(p, (int)len, &filtered) < 0 ||
          kputc('\n', &filtered) < 0){
        free(filtered.s);
        return NULL;
      }
    }
    if (!e) break;
    p = e + 1;
  }
  out = sam_hdr_parse(filtered.l, filtered.s ? filtered.s : "");
  free(filtered.s);
  return out;
}

static int merge_header_line(bam_hdr_t* dst, const char* line, size_t len,
                             const pg_id_map_t* pg_map, int strip_pg){
  kstring_t existing = {0, 0, NULL};
  kstring_t rewritten = {0, 0, NULL};
  char* id = NULL;
  char type[3];
  int find_rc;
  int ok = 0;

  if (len == 0 || line[0] != '@') return 0;
  if (line_type_is(line, len, "HD") || line_type_is(line, len, "SQ")) return 0;
  if (strip_pg && line_type_is(line, len, "PG")) return 0;

  if (line_type_is(line, len, "PG") && pg_map){
    if (rewrite_pg_header_line(line, len, pg_map, &rewritten) != 0){
      free(rewritten.s);
      return -1;
    }
    ok = merge_header_line(dst, rewritten.s, rewritten.l, NULL, 0);
    free(rewritten.s);
    return ok;
  }

  if (line_type_is(line, len, "RG") || line_type_is(line, len, "PG")){
    type[0] = line[1];
    type[1] = line[2];
    type[2] = '\0';
    id = extract_header_tag_value_dup(line, len, "ID");
    if (!id){
      LOG_PIPELINE("header line missing ID: %.*s", (int)len, line);
      return -1;
    }
    find_rc = sam_hdr_find_line_id(dst, type, "ID", id, &existing);
    if (find_rc == 0){
      if (existing.l == len && strncmp(existing.s, line, len) == 0){
        ok = 0;
      } else {
        LOG_PIPELINE("conflicting @%c%c header line for ID=%s", type[0], type[1], id);
        ok = -1;
      }
    } else if (find_rc == -1){
      ok = sam_hdr_add_lines(dst, line, len);
      if (ok != 0) LOG_PIPELINE("failed to append @%c%c header line for ID=%s", type[0], type[1], id);
    } else {
      LOG_PIPELINE("failed to inspect @%c%c header line for ID=%s", type[0], type[1], id);
      ok = -1;
    }
    free(id);
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

static int merge_header_metadata(bam_hdr_t* dst, const bam_hdr_t* src,
                                 const pg_id_map_t* pg_map, int strip_pg){
  const char* text = sam_hdr_str((bam_hdr_t*)src);
  const char* p = text;
  if (!text) return -1;
  while (*p){
    const char* e = strchr(p, '\n');
    size_t len = e ? (size_t)(e - p) : strlen(p);
    if (merge_header_line(dst, p, len, pg_map, strip_pg) != 0) return -1;
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int build_merged_header(const cli_opts_t* o, bam_hdr_t** out_hdr,
                               pg_id_map_t** out_pg_maps){
  io_ctx_t io0 = {0};
  bam_hdr_t* merged = NULL;
  pg_id_map_t* pg_maps = (pg_id_map_t*)calloc(o->bam_list.n, sizeof(*pg_maps));
  strvec_t used_pg_ids;
  strvec_init(&used_pg_ids);
  if (!pg_maps) return -1;
  if (io_open(o->bam_list.data[0], NULL, &io0) != 0){ free(pg_maps); return -1; }
  merged = o->strip_pg ? duplicate_header_without_pg(io0.hdr) : sam_hdr_dup(io0.hdr);
  if (!merged || (!o->strip_pg &&
                  (build_pg_id_map(io0.hdr, 0, &used_pg_ids, &pg_maps[0]) != 0 ||
                   validate_pg_parents(io0.hdr, 0, &pg_maps[0]) != 0))){
    bam_hdr_destroy(merged);
    io_close(&io0);
    strvec_free(&used_pg_ids);
    free_pg_id_maps(pg_maps, o->bam_list.n);
    return -1;
  }
  io_close(&io0);

  for (size_t fi = 1; fi < o->bam_list.n; ++fi){
    io_ctx_t io = {0};
    if (io_open(o->bam_list.data[fi], NULL, &io) != 0){
      bam_hdr_destroy(merged); strvec_free(&used_pg_ids); free_pg_id_maps(pg_maps, o->bam_list.n); return -1;
    }
    if (!headers_compatible(merged, io.hdr)){
      LOG_PIPELINE("[split] header mismatch in %s", o->bam_list.data[fi]);
      io_close(&io);
      bam_hdr_destroy(merged);
      strvec_free(&used_pg_ids);
      free_pg_id_maps(pg_maps, o->bam_list.n);
      return -1;
    }
    if ((!o->strip_pg &&
         (build_pg_id_map(io.hdr, fi, &used_pg_ids, &pg_maps[fi]) != 0 ||
          validate_pg_parents(io.hdr, fi, &pg_maps[fi]) != 0)) ||
        adopt_optional_sq_tags(merged, io.hdr) != 0 ||
        merge_header_metadata(merged, io.hdr,
                              o->strip_pg ? NULL : &pg_maps[fi], o->strip_pg) != 0){
      io_close(&io);
      bam_hdr_destroy(merged);
      strvec_free(&used_pg_ids);
      free_pg_id_maps(pg_maps, o->bam_list.n);
      return -1;
    }
    io_close(&io);
  }

  if ((!o->strip_pg && add_isoumi_program_record(merged, o, &used_pg_ids) != 0) ||
      mark_header_sort_unknown(merged) != 0){
    bam_hdr_destroy(merged);
    strvec_free(&used_pg_ids);
    free_pg_id_maps(pg_maps, o->bam_list.n);
    return -1;
  }
  strvec_free(&used_pg_ids);
  *out_hdr = merged;
  *out_pg_maps = pg_maps;
  return 0;
}

static void cleanup_tmp(const cli_opts_t* o);

static int phase_split(const cli_opts_t* o){
  bam_hdr_t* merged_hdr = NULL;
  pg_id_map_t* pg_maps = NULL;
  samFile** outs = NULL;
  bam1_t* b = NULL;
  int opened_writers = 0;
  int result = -1;
  char path[4096];
  if (ensure_dir(o->tmp_dir)!=0) return -1;
  if (build_merged_header(o, &merged_hdr, &pg_maps) != 0){
    if (!o->keep_tmp) cleanup_tmp(o);
    return -1;
  }
  outs = (samFile**)calloc((size_t)o->buckets, sizeof(samFile*));
  if (!outs) goto cleanup;
  for (int i=0;i<o->buckets;++i){
    snprintf(path, sizeof(path), "%s/bucket_%05d.bam", o->tmp_dir, i);
    /* Bucket BAMs are transient and will be recompressed into the final BAM. */
    outs[i] = sam_open(path, "wb1");
    if (!outs[i]){ LOG_PIPELINE("[split] cannot create %s", path); goto cleanup; }
    opened_writers = i + 1;
    if (sam_hdr_write(outs[i], merged_hdr) < 0){ LOG_PIPELINE("[split] hdr write fail %s", path); goto cleanup; }
  }
  b = bam_init1();
  if (!b) goto cleanup;
  for (size_t fi=0; fi<o->bam_list.n; ++fi){
    const char* inbam = o->bam_list.data[fi];
    char input_scope[64];
    samFile* in = NULL;
    bam_hdr_t* hdr = NULL;
    int read_rc = 0;
    int failed = 0;
    if (input_scope_label(fi, input_scope, sizeof(input_scope)) != 0){
      LOG_PIPELINE("[split] failed to build input scope label");
      goto cleanup;
    }
    in = sam_open(inbam, "r");
    if (!in){ LOG_PIPELINE("[split] cannot open %s", inbam); goto cleanup; }
    if (o->threads > 1 && hts_set_threads(in, o->threads) < 0){
      LOG_PIPELINE("[split] warning: cannot enable threaded input decoding for %s", inbam);
    }
    hdr = sam_hdr_read(in);
    if (!hdr){ LOG_PIPELINE("[split] cannot read header %s", inbam); sam_close(in); goto cleanup; }
    if (!headers_compatible(merged_hdr, hdr)){
      LOG_PIPELINE("[split] header mismatch in %s", inbam);
      bam_hdr_destroy(hdr); sam_close(in); goto cleanup;
    }
    while ((read_rc = sam_read1(in, hdr, b)) >= 0){
      if ((o->strip_pg ? strip_record_pg(b, inbam)
                       : rewrite_record_pg(b, &pg_maps[fi], inbam)) != 0 ||
          set_input_scope_tag(o, b, input_scope) != 0){ failed = 1; break; }
      const char* cb = get_tag_Z(b, o->cell_tag);
      if (is_non_primary(b) && !is_unmapped(b) && !cb){
        LOG_PIPELINE("[split] mapped non-primary record %s in %s lacks %s; "
                     "cannot reliably inherit its primary assignment",
                     bam_get_qname(b), inbam, o->cell_tag);
        failed = 1;
        break;
      }
      unsigned long bid = hash_cb(cb) % (unsigned long)o->buckets;
      if (sam_write1(outs[bid], merged_hdr, b) < 0){ LOG_PIPELINE("[split] write fail for bucket %lu", bid); failed = 1; break; }
    }
    if (read_rc < -1){ LOG_PIPELINE("[split] read fail %s", inbam); failed = 1; }
    bam_hdr_destroy(hdr); sam_close(in);
    if (failed) goto cleanup;
  }
  result = 0;

cleanup:
  bam_destroy1(b);
  bam_hdr_destroy(merged_hdr);
  if (close_writers(outs, opened_writers) != 0){
    LOG_PIPELINE("[split] failed to flush/close one or more bucket files");
    result = -1;
  }
  free(outs);
  free_pg_id_maps(pg_maps, o->bam_list.n);
  if (result != 0 && !o->keep_tmp) cleanup_tmp(o);
  return result;
}

typedef struct {
  char* key;
  char* umi;
  int count;
  double* pos_qsum;
  int quality_count;
  int umi_len;
  int has_qual;
  long best_ordinal;
  int best_mapq;
  int best_ref_span;
  int best_query_len;
  uint64_t hash;
} aggregate_t;

typedef struct {
  aggregate_t* slots;
  size_t n;
  size_t cap;
} aggregate_table_t;

static int cmp_aggregate(const void* a, const void* b){
  const aggregate_t* x=(const aggregate_t*)a;
  const aggregate_t* y=(const aggregate_t*)b;
  int c=strcmp(x->key, y->key);
  if(c) return c;
  return strcmp(x->umi, y->umi);
}

static int read_ref_span(const bam1_t* b);
static int build_read_identity_buf(const cli_opts_t* o, bam1_t* b,
                                   kstring_t* out);
typedef struct {
  char* umi;
  int count;
  double mean_qual;
  double* pos_qsum;
  int quality_count;
  int umi_len;
  int has_qual;
} umi_stat_t;
typedef struct {
  int rep_idx;
  int parent_idx;
  int hamming;
  int edge_hamming;
  int path_length;
  double corr_mismatch_q;
  double raw_mismatch_q;
  double confidence;
  int mismatch_q_n;
} umi_assignment_t;
static int cmp_umi_stat_desc(const void* a, const void* b){
  const umi_stat_t* x=(const umi_stat_t*)a;
  const umi_stat_t* y=(const umi_stat_t*)b;
  if (x->count!=y->count) return (y->count-x->count);
  return strcmp(x->umi, y->umi);
}
typedef struct {
  char* raw_umi;
  char* corr_umi;
  char* parent_umi;
  int raw_count;
  int seed_count;
  int parent_count;
  int hamming;
  int edge_hamming;
  int path_length;
  double raw_avgq;
  double seed_avgq;
  double confidence;
  int quality_supported;
  correction_method_t method;
} map_item_t;
static int cmp_map_item(const void* a, const void* b){ const map_item_t* x=(const map_item_t*)a; const map_item_t* y=(const map_item_t*)b; int c=strcmp(x->raw_umi, y->raw_umi); if(c) return c; return strcmp(x->corr_umi, y->corr_umi); }
typedef struct { char* key; map_item_t* items; int n; } key_map_t;
static int cmp_key_map(const void* a, const void* b){ const key_map_t* x=(const key_map_t*)a; const key_map_t* y=(const key_map_t*)b; return strcmp(x->key, y->key); }
static char* sdup(const char* s){ size_t n=s?strlen(s):0; char* p=(char*)malloc(n+1); if(!p) return NULL; memcpy(p,s?s:"",n+1); return p; }
static void free_umi_stats(umi_stat_t* stats, int n){
  for (int i=0;i<n;++i){
    free(stats[i].pos_qsum);
  }
  free(stats);
}
static void free_aggregates(aggregate_t* arr, size_t n){
  if (!arr) return;
  for (size_t i=0; i<n; ++i){
    free(arr[i].key);
    free(arr[i].umi);
    free(arr[i].pos_qsum);
  }
  free(arr);
}
static void free_map_items(map_item_t* items, int n){ for (int i=0;i<n;++i){ free(items[i].raw_umi); free(items[i].corr_umi); free(items[i].parent_umi); } free(items); }
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

static void finalize_umi_quality(umi_stat_t* st){
  int i;
  double sum = 0.0;
  int64_t n = 0;
  if (!st->has_qual || !st->pos_qsum || st->quality_count <= 0){
    st->mean_qual = -1.0;
    return;
  }
  for (i=0; i<st->umi_len; ++i){
    sum += st->pos_qsum[i];
    n += st->quality_count;
  }
  st->mean_qual = n > 0 ? sum / (double)n : -1.0;
}

static uint64_t hash_key_umi(const char* key, const char* umi){
  uint64_t h = UINT64_C(14695981039346656037);
  const unsigned char* p;
  for (p=(const unsigned char*)key; *p; ++p){ h ^= *p; h *= UINT64_C(1099511628211); }
  h ^= UINT64_C(0xff);
  h *= UINT64_C(1099511628211);
  for (p=(const unsigned char*)umi; *p; ++p){ h ^= *p; h *= UINT64_C(1099511628211); }
  return h ? h : UINT64_C(1);
}

static int aggregate_table_resize(aggregate_table_t* table, size_t new_cap){
  aggregate_t* grown;
  if (!table || new_cap < 16 || (new_cap & (new_cap - 1)) != 0) return -1;
  grown = (aggregate_t*)calloc(new_cap, sizeof(*grown));
  if (!grown) return -1;
  for (size_t i=0; i<table->cap; ++i){
    aggregate_t* old = &table->slots[i];
    size_t idx;
    if (!old->key) continue;
    idx = (size_t)old->hash & (new_cap - 1);
    while (grown[idx].key) idx = (idx + 1) & (new_cap - 1);
    grown[idx] = *old;
  }
  free(table->slots);
  table->slots = grown;
  table->cap = new_cap;
  return 0;
}

static int aggregate_add_quality(aggregate_t* entry, const char* qual){
  if (!qual || !entry || !entry->umi) return 0;
  if ((int)strlen(qual) != entry->umi_len) return 0;
  if (!entry->pos_qsum){
    entry->pos_qsum = (double*)calloc((size_t)entry->umi_len, sizeof(double));
    if (!entry->pos_qsum){
      free(entry->pos_qsum);
      entry->pos_qsum = NULL;
      return -1;
    }
  }
  for (int i=0; i<entry->umi_len; ++i){
    entry->pos_qsum[i] += (double)((int)(unsigned char)qual[i] - 33);
  }
  entry->quality_count++;
  entry->has_qual = 1;
  return 0;
}

static int summary_is_better(int mapq, int ref_span, int query_len, long ordinal,
                             const aggregate_t* current){
  if (current->best_ordinal < 0) return 1;
  if (mapq != current->best_mapq) return mapq > current->best_mapq;
  if (ref_span != current->best_ref_span) return ref_span > current->best_ref_span;
  if (query_len != current->best_query_len) return query_len > current->best_query_len;
  return ordinal < current->best_ordinal;
}

static void aggregate_set_best(aggregate_t* entry, const bam1_t* b, long ordinal){
  entry->best_ordinal = ordinal;
  entry->best_mapq = b ? b->core.qual : 0;
  entry->best_ref_span = read_ref_span(b);
  entry->best_query_len = b ? b->core.l_qseq : 0;
}

static int aggregate_table_add(aggregate_table_t* table, char* key,
                               const char* umi, const char* qual,
                               const bam1_t* b, long ordinal){
  uint64_t hash;
  size_t idx;
  aggregate_t* entry;
  if (!table || !key || !umi){ free(key); return -1; }
  if (table->cap == 0 && aggregate_table_resize(table, 1024) != 0){ free(key); return -1; }
  if (table->n + 1 >= table->cap - table->cap / 4){
    if (table->cap > SIZE_MAX / 2 ||
        aggregate_table_resize(table, table->cap * 2) != 0){ free(key); return -1; }
  }
  hash = hash_key_umi(key, umi);
  idx = (size_t)hash & (table->cap - 1);
  while (table->slots[idx].key){
    entry = &table->slots[idx];
    if (entry->hash == hash && strcmp(entry->key, key) == 0 &&
        strcmp(entry->umi, umi) == 0){
      int mapq = b ? b->core.qual : 0;
      int ref_span = read_ref_span(b);
      int query_len = b ? b->core.l_qseq : 0;
      free(key);
      if (entry->count == INT_MAX || aggregate_add_quality(entry, qual) != 0) return -1;
      entry->count++;
      if (summary_is_better(mapq, ref_span, query_len, ordinal, entry)){
        entry->best_ordinal = ordinal;
        entry->best_mapq = mapq;
        entry->best_ref_span = ref_span;
        entry->best_query_len = query_len;
      }
      return 0;
    }
    idx = (idx + 1) & (table->cap - 1);
  }
  entry = &table->slots[idx];
  entry->key = key;
  entry->umi = sdup(umi);
  if (!entry->umi){
    free(entry->key);
    entry->key = NULL;
    return -1;
  }
  entry->count = 1;
  entry->umi_len = (int)strlen(umi);
  entry->best_ordinal = -1;
  entry->hash = hash;
  aggregate_set_best(entry, b, ordinal);
  if (aggregate_add_quality(entry, qual) != 0){
    free(entry->key);
    free(entry->umi);
    memset(entry, 0, sizeof(*entry));
    return -1;
  }
  table->n++;
  return 0;
}

static void aggregate_table_destroy(aggregate_table_t* table){
  if (!table) return;
  if (table->slots){
    for (size_t i=0; i<table->cap; ++i){
      aggregate_t* entry = &table->slots[i];
      if (!entry->key) continue;
      free(entry->key);
      free(entry->umi);
      free(entry->pos_qsum);
    }
  }
  free(table->slots);
  table->slots = NULL;
  table->n = table->cap = 0;
}

static aggregate_t* aggregate_table_take_sorted(aggregate_table_t* table){
  aggregate_t* dense;
  size_t out = 0;
  if (!table || table->n == 0){
    aggregate_table_destroy(table);
    return NULL;
  }
  dense = (aggregate_t*)malloc(sizeof(*dense) * table->n);
  if (!dense) return NULL;
  for (size_t i=0; i<table->cap; ++i){
    if (table->slots[i].key) dense[out++] = table->slots[i];
  }
  free(table->slots);
  table->slots = NULL;
  table->cap = table->n = 0;
  qsort(dense, out, sizeof(*dense), cmp_aggregate);
  return dense;
}

static double pos_mean_qual(const umi_stat_t* st, int pos){
  if (!st->has_qual || !st->pos_qsum || st->quality_count <= 0) return -1.0;
  if (pos < 0 || pos >= st->umi_len) return -1.0;
  return st->pos_qsum[pos] / (double)st->quality_count;
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

static int can_assign_to_seed(const cli_opts_t* o, const umi_stat_t* seed, const umi_stat_t* raw,
                              umi_assignment_t* out){
  int ham = 0, mismatch_q_n = 0;
  double seed_mq = -1.0, raw_mq = -1.0;
  double conf;

  if (!seed || !raw || seed == raw) return 0;
  if (seed->count <= 0 || seed->count < raw->count) return 0;
  if (seed->count == raw->count && !better_representative(o, seed, raw)) return 0;
  if ((double)raw->count / (double)seed->count > o->ratio) return 0;
  if (compute_umi_distance_metrics(seed, raw, o->ham, &ham, &seed_mq, &raw_mq, &mismatch_q_n) != 0) return 0;

  conf = compute_merge_confidence(o, seed, raw, ham, seed_mq, raw_mq, mismatch_q_n);
  if (o->min_merge_confidence > 0.0 && conf < o->min_merge_confidence) return 0;

  if (out){
    out->hamming = ham;
    out->corr_mismatch_q = seed_mq;
    out->raw_mismatch_q = raw_mq;
    out->confidence = conf;
    out->mismatch_q_n = mismatch_q_n;
  }
  return 1;
}

static int choose_next_seed(const cli_opts_t* o, const umi_stat_t* stats, int n,
                            const umi_assignment_t* assignments){
  int best = -1;
  for (int i=0; i<n; ++i){
    if (assignments[i].rep_idx >= 0) continue;
    if (best < 0 || better_representative(o, &stats[i], &stats[best])) best = i;
  }
  return best;
}

static int build_ratio_assignments(const cli_opts_t* o, const umi_stat_t* stats, int n,
                                   umi_assignment_t** out){
  umi_assignment_t* assignments = (umi_assignment_t*)calloc((size_t)n, sizeof(umi_assignment_t));
  if (!assignments) return -1;
  for (int i=0; i<n; ++i){
    assignments[i].rep_idx = -1;
    assignments[i].parent_idx = -1;
  }

  for (;;) {
    int seed = choose_next_seed(o, stats, n, assignments);
    if (seed < 0) break;

    assignments[seed].rep_idx = seed;
    assignments[seed].parent_idx = seed;
    assignments[seed].hamming = 0;
    assignments[seed].edge_hamming = 0;
    assignments[seed].path_length = 0;
    assignments[seed].corr_mismatch_q = -1.0;
    assignments[seed].raw_mismatch_q = -1.0;
    assignments[seed].confidence = 1.0;
    assignments[seed].mismatch_q_n = 0;

    for (int idx=0; idx<n; ++idx){
      umi_assignment_t candidate = {0};
      if (assignments[idx].rep_idx >= 0 || idx == seed) continue;
      if (can_assign_to_seed(o, &stats[seed], &stats[idx], &candidate)){
        candidate.rep_idx = seed;
        candidate.parent_idx = seed;
        candidate.edge_hamming = candidate.hamming;
        candidate.path_length = 1;
        assignments[idx] = candidate;
      }
    }
  }

  *out = assignments;
  return 0;
}

typedef struct {
  int* data;
  int n;
  int cap;
} int_vec_t;

static int int_vec_push(int_vec_t* v, int value){
  if (v->n == v->cap){
    int new_cap = v->cap ? v->cap * 2 : 4;
    int* new_data = (int*)realloc(v->data, sizeof(int) * (size_t)new_cap);
    if (!new_data) return -1;
    v->data = new_data;
    v->cap = new_cap;
  }
  v->data[v->n++] = value;
  return 0;
}

static void free_directional_graph(int_vec_t* graph, int n){
  if (!graph) return;
  for (int i=0; i<n; ++i) free(graph[i].data);
  free(graph);
}

static int directional_count_allows(const umi_stat_t* from, const umi_stat_t* to){
  return (int64_t)from->count >= (int64_t)2 * (int64_t)to->count - 1;
}

static int build_directional_graph(const cli_opts_t* o, const umi_stat_t* stats, int n,
                                   int_vec_t** out){
  int_vec_t* graph = (int_vec_t*)calloc((size_t)n, sizeof(int_vec_t));
  if (!graph) return -1;

  for (int i=0; i<n; ++i){
    for (int j=i+1; j<n; ++j){
      int ham = 0, mismatch_q_n = 0;
      double i_mq = -1.0, j_mq = -1.0;
      if (compute_umi_distance_metrics(&stats[i], &stats[j], o->ham,
                                       &ham, &i_mq, &j_mq, &mismatch_q_n) != 0){
        continue;
      }
      if (directional_count_allows(&stats[i], &stats[j]) &&
          int_vec_push(&graph[i], j) != 0){
        free_directional_graph(graph, n);
        return -1;
      }
      if (directional_count_allows(&stats[j], &stats[i]) &&
          int_vec_push(&graph[j], i) != 0){
        free_directional_graph(graph, n);
        return -1;
      }
    }
  }

  *out = graph;
  return 0;
}

static int build_directional_assignments(const cli_opts_t* o, const umi_stat_t* stats, int n,
                                         umi_assignment_t** out){
  int_vec_t* graph = NULL;
  int* queue = NULL;
  umi_assignment_t* assignments = NULL;

  if (build_directional_graph(o, stats, n, &graph) != 0) return -1;
  assignments = (umi_assignment_t*)calloc((size_t)n, sizeof(umi_assignment_t));
  queue = (int*)malloc(sizeof(int) * (size_t)n);
  if (!assignments || !queue){
    free(assignments);
    free(queue);
    free_directional_graph(graph, n);
    return -1;
  }
  for (int i=0; i<n; ++i){
    assignments[i].rep_idx = -1;
    assignments[i].parent_idx = -1;
  }

  for (int seed=0; seed<n; ++seed){
    int head = 0, tail = 0;
    if (assignments[seed].rep_idx >= 0) continue;

    assignments[seed].rep_idx = seed;
    assignments[seed].parent_idx = seed;
    assignments[seed].hamming = 0;
    assignments[seed].edge_hamming = 0;
    assignments[seed].path_length = 0;
    assignments[seed].corr_mismatch_q = -1.0;
    assignments[seed].raw_mismatch_q = -1.0;
    assignments[seed].confidence = 1.0;
    assignments[seed].mismatch_q_n = 0;
    queue[tail++] = seed;

    while (head < tail){
      int parent = queue[head++];
      for (int k=0; k<graph[parent].n; ++k){
        int child = graph[parent].data[k];
        int edge_ham = 0, edge_q_n = 0;
        int root_ham = 0, root_q_n = 0;
        double parent_mq = -1.0, child_mq = -1.0;
        double root_mq = -1.0, root_child_mq = -1.0;
        if (assignments[child].rep_idx >= 0) continue;
        if (compute_umi_distance_metrics(&stats[parent], &stats[child], o->ham,
                                         &edge_ham, &parent_mq, &child_mq,
                                         &edge_q_n) != 0 ||
            compute_umi_distance_metrics(&stats[seed], &stats[child], -1,
                                         &root_ham, &root_mq, &root_child_mq,
                                         &root_q_n) != 0){
          free(assignments);
          free(queue);
          free_directional_graph(graph, n);
          return -1;
        }
        assignments[child].rep_idx = seed;
        assignments[child].parent_idx = parent;
        assignments[child].hamming = root_ham;
        assignments[child].edge_hamming = edge_ham;
        assignments[child].path_length = assignments[parent].path_length + 1;
        assignments[child].corr_mismatch_q = parent_mq;
        assignments[child].raw_mismatch_q = child_mq;
        assignments[child].confidence = compute_merge_confidence(
            o, &stats[parent], &stats[child], edge_ham,
            parent_mq, child_mq, edge_q_n);
        assignments[child].mismatch_q_n = edge_q_n;
        queue[tail++] = child;
      }
    }
  }

  free(queue);
  free_directional_graph(graph, n);
  *out = assignments;
  return 0;
}

static int build_umi_assignments(const cli_opts_t* o, const umi_stat_t* stats, int n,
                                 umi_assignment_t** out){
  if (o->correction_method == CORRECTION_DIRECTIONAL){
    return build_directional_assignments(o, stats, n, out);
  }
  return build_ratio_assignments(o, stats, n, out);
}

static const char* correction_method_name(correction_method_t method){
  return method == CORRECTION_DIRECTIONAL ? "directional" : "ratio";
}

static const char* merge_reason(const map_item_t* item){
  if (!item) return "unknown";
  if (item->hamming == 0 || strcmp(item->raw_umi, item->corr_umi) == 0) return "self";
  if (item->method == CORRECTION_DIRECTIONAL) return "directional";
  if (item->quality_supported) return "count+quality";
  return "count+distance";
}

static int emit_correction_row(FILE* fp, const char* key, const map_item_t* item, int bucket){
  if (!fp || !item) return 0;
  if (fprintf(fp, "%s\t%s\t%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.4f\t%s\t%d\t%s\t%s\t%d\t%d\t%d\n",
              key, item->raw_umi, item->corr_umi, item->raw_count, item->seed_count,
              item->hamming, item->raw_avgq, item->seed_avgq, item->confidence,
              merge_reason(item), bucket, correction_method_name(item->method),
              item->parent_umi, item->parent_count, item->edge_hamming,
              item->path_length) < 0){
    return -1;
  }
  return 0;
}

static int build_bucket_mapping(const cli_opts_t* o, const char* bucket_bam,
                                FILE* explain_fp, int bucket_id,
                                key_map_t** out_maps, int* out_nm,
                                aggregate_t** out_aggregates, size_t* out_n_aggregates,
                                strvec_t* non_primary_keys){
  io_ctx_t io = {0};
  bam1_t* b = NULL;
  aggregate_table_t table = {0};
  aggregate_t* arr = NULL;
  key_map_t* km = NULL;
  int km_n = 0, km_cap = 64;
  int use_qual = want_umi_quality(o);
  int read_rc = 0;
  long ordinal = 0;
  unsigned long long primary_reads = 0;
  size_t n = 0;
  kstring_t read_identity = {0, 0, NULL};

  string_set_t non_primary_set = {0};
  if (io_open(bucket_bam, NULL, &io) != 0) return -1;
  b = bam_init1();
  if (!b) goto fail;
  while ((read_rc = sam_read1(io.in, io.hdr, b)) >= 0){
    long current_ordinal = ordinal++;
    const char* umi;
    const char* cb;
    const char* umi_qual = NULL;
    char* key;
    if (is_unmapped(b)) continue;
    if (is_non_primary(b)){
      if (build_read_identity_buf(o, b, &read_identity) != 0 ||
          string_set_add(&non_primary_set, read_identity.s) != 0) goto fail;
      continue;
    }
    umi = get_tag_Z(b, o->umi_tag);
    cb = get_tag_Z(b, o->cell_tag);
    if (!umi || !cb) continue;
    key = build_group_key(b, o->cell_tag, o->gene_tag,
                          o->source_tag, active_input_scope_tag(o),
                          o->no_gene, o->no_structure,
                          o->locus_bin, o->sj_jitter, o->end_bin);
    if (!key) goto fail;
    if (use_qual) umi_qual = get_tag_Z(b, o->umi_qual_tag);
    if (aggregate_table_add(&table, key, umi, umi_qual, b, current_ordinal) != 0){
      goto fail;
    }
    primary_reads++;
  }
  if (read_rc < -1) goto fail;
  bam_destroy1(b);
  b = NULL;
  io_close(&io);
  n = table.n;
  arr = aggregate_table_take_sorted(&table);
  if (n > 0 && !arr) goto fail;
  if (primary_reads > 0){
    LOG_BUCKET(bucket_id, "aggregated %llu primary records into %zu unique group/UMI entries",
               primary_reads, n);
  }

  km = (key_map_t*)malloc(sizeof(*km) * (size_t)km_cap);
  if (!km) goto fail;
  for (size_t i=0; i<n; ){
    size_t j = i + 1;
    umi_stat_t* uc;
    umi_assignment_t* assignments = NULL;
    map_item_t* items;
    int uc_n;
    int mi_n = 0;
    while (j<n && strcmp(arr[j].key, arr[i].key)==0) j++;
    if (j - i > (size_t)INT_MAX) goto fail;
    uc_n = (int)(j - i);
    uc = (umi_stat_t*)calloc((size_t)uc_n, sizeof(*uc));
    if (!uc) goto fail;
    for (int idx=0; idx<uc_n; ++idx){
      aggregate_t* entry = &arr[i + (size_t)idx];
      uc[idx].umi = entry->umi;
      uc[idx].count = entry->count;
      uc[idx].umi_len = entry->umi_len;
      uc[idx].has_qual = entry->has_qual;
      uc[idx].quality_count = entry->quality_count;
      uc[idx].pos_qsum = entry->pos_qsum;
      entry->pos_qsum = NULL;
      finalize_umi_quality(&uc[idx]);
    }
    qsort(uc, (size_t)uc_n, sizeof(*uc), cmp_umi_stat_desc);
    if (build_umi_assignments(o, uc, uc_n, &assignments) != 0){
      free_umi_stats(uc, uc_n);
      goto fail;
    }
    items = (map_item_t*)calloc((size_t)uc_n, sizeof(*items));
    if (!items){
      free(assignments);
      free_umi_stats(uc, uc_n);
      goto fail;
    }
    for (int idx=0; idx<uc_n; ++idx){
      int ridx=assignments[idx].rep_idx;
      int pidx=assignments[idx].parent_idx;
      items[mi_n].raw_umi=sdup(uc[idx].umi);
      items[mi_n].corr_umi=sdup(uc[ridx].umi);
      items[mi_n].parent_umi=sdup(uc[pidx].umi);
      if (!items[mi_n].raw_umi || !items[mi_n].corr_umi || !items[mi_n].parent_umi){
        free_map_items(items, mi_n + 1);
        free(assignments);
        free_umi_stats(uc, uc_n);
        goto fail;
      }
      items[mi_n].raw_count = uc[idx].count;
      items[mi_n].seed_count = uc[ridx].count;
      items[mi_n].parent_count = uc[pidx].count;
      items[mi_n].raw_avgq = uc[idx].mean_qual;
      items[mi_n].seed_avgq = uc[ridx].mean_qual;
      items[mi_n].hamming = assignments[idx].hamming;
      items[mi_n].edge_hamming = assignments[idx].edge_hamming;
      items[mi_n].path_length = assignments[idx].path_length;
      items[mi_n].method = o->correction_method;
      items[mi_n].quality_supported = o->correction_method == CORRECTION_RATIO &&
        o->quality_aware && assignments[idx].mismatch_q_n > 0 &&
        assignments[idx].corr_mismatch_q >= 0.0 &&
        assignments[idx].raw_mismatch_q >= 0.0 &&
        assignments[idx].corr_mismatch_q > assignments[idx].raw_mismatch_q;
      items[mi_n].confidence = assignments[idx].confidence;
      if (emit_correction_row(explain_fp, arr[i].key, &items[mi_n], bucket_id) != 0){
        free_map_items(items, mi_n + 1);
        free(assignments);
        free_umi_stats(uc, uc_n);
        goto fail;
      }
      mi_n++;
    }
    qsort(items, (size_t)mi_n, sizeof(*items), cmp_map_item);
    if (km_n == km_cap){
      int nc;
      key_map_t* grown;
      if (km_cap > INT_MAX / 2){
        free_map_items(items, mi_n);
        free(assignments);
        free_umi_stats(uc, uc_n);
        goto fail;
      }
      nc = km_cap * 2;
      grown = (key_map_t*)realloc(km, sizeof(*grown) * (size_t)nc);
      if (!grown){
        free_map_items(items, mi_n);
        free(assignments);
        free_umi_stats(uc, uc_n);
        goto fail;
      }
      km = grown;
      km_cap = nc;
    }
    km[km_n].key=sdup(arr[i].key);
    if (!km[km_n].key){
      free_map_items(items, mi_n);
      free(assignments);
      free_umi_stats(uc, uc_n);
      goto fail;
    }
    km[km_n].items=items;
    km[km_n].n=mi_n;
    km_n++;
    free(assignments);
    free_umi_stats(uc, uc_n);
    i=j;
  }
  qsort(km, (size_t)km_n, sizeof(*km), cmp_key_map);
  if (string_set_export(&non_primary_set, non_primary_keys) != 0) goto fail;
  string_set_destroy(&non_primary_set);
  *out_maps = km;
  *out_nm = km_n;
  *out_aggregates = arr;
  *out_n_aggregates = n;
  free(read_identity.s);
  return 0;

fail:
  free(read_identity.s);
  bam_destroy1(b);
  string_set_destroy(&non_primary_set);
  io_close(&io);
  aggregate_table_destroy(&table);
  free_aggregates(arr, n);
  free_key_maps(km, km_n);
  return -1;
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

typedef struct {
  char* key;
  char* corr;
  int count;
  long best_ordinal;
  int best_mapq;
  int best_ref_span;
  int best_query_len;
} molecule_t;
static int molecule_cmp_key(const molecule_t* x, const molecule_t* y){
  int c=strcmp(x->key,y->key);
  if(c) return c;
  return strcmp(x->corr,y->corr);
}
static int molecule_find(molecule_t* arr, int n, const char* key, const char* corr){
  int lo=0, hi=n-1;
  molecule_t tmp={(char*)key,(char*)corr,0,0,0,0,0};
  while(lo<=hi){
    int mid=(lo+hi)/2;
    int c=molecule_cmp_key(&tmp,&arr[mid]);
    if(c==0) return mid;
    if(c<0) hi=mid-1; else lo=mid+1;
  }
  return -lo-1;
}
static void free_molecules(molecule_t* mols, int mol_n){
  for (int s=0;s<mol_n;++s){ free(mols[s].key); free(mols[s].corr); }
  free(mols);
}
static char* build_molecule_id(const char* cb, const char* key, const char* corr){
  char* cb_escaped = cb ? escape_key_component(cb) : NULL;
  char* corr_escaped = corr ? escape_key_component(corr) : NULL;
  const char* cbv = cb ? cb_escaped : "NA";
  const char* keyv = key ? key : "NA";
  const char* corrv = corr ? corr_escaped : "NA";
  char* out;
  int need;
  if ((cb && !cb_escaped) || (corr && !corr_escaped)){
    free(cb_escaped);
    free(corr_escaped);
    return NULL;
  }
  need = snprintf(NULL, 0, "%s|%s|%s", cbv, keyv, corrv);
  if (need < 0){
    free(cb_escaped); free(corr_escaped);
    return NULL;
  }
  out = (char*)malloc((size_t)need + 1);
  if (!out){
    free(cb_escaped); free(corr_escaped);
    return NULL;
  }
  snprintf(out, (size_t)need + 1, "%s|%s|%s", cbv, keyv, corrv);
  free(cb_escaped); free(corr_escaped);
  return out;
}

static int read_ref_span(const bam1_t* b){
  hts_pos_t end_pos;
  hts_pos_t span;
  if (!b || b->core.pos < 0) return 0;
  end_pos = bam_endpos(b);
  span = end_pos > b->core.pos ? end_pos - b->core.pos : 0;
  return span > INT_MAX ? INT_MAX : (int)span;
}

typedef struct {
  const char* key;
  const char* corr;
  int count;
  long best_ordinal;
  int best_mapq;
  int best_ref_span;
  int best_query_len;
} molecule_candidate_t;

static int cmp_molecule_candidate(const void* a, const void* b){
  const molecule_candidate_t* x = (const molecule_candidate_t*)a;
  const molecule_candidate_t* y = (const molecule_candidate_t*)b;
  int c = strcmp(x->key, y->key);
  if (c) return c;
  return strcmp(x->corr, y->corr);
}

static int candidate_is_better(const molecule_candidate_t* candidate,
                               const molecule_t* current){
  if (current->best_ordinal < 0) return 1;
  if (candidate->best_mapq != current->best_mapq)
    return candidate->best_mapq > current->best_mapq;
  if (candidate->best_ref_span != current->best_ref_span)
    return candidate->best_ref_span > current->best_ref_span;
  if (candidate->best_query_len != current->best_query_len)
    return candidate->best_query_len > current->best_query_len;
  return candidate->best_ordinal < current->best_ordinal;
}

static int build_bucket_molecules_from_aggregates(const key_map_t* km, int km_n,
                                                   const aggregate_t* aggregates,
                                                   size_t n_aggregates,
                                                   molecule_t** out_mols, int* out_n){
  molecule_candidate_t* candidates = NULL;
  molecule_t* mols = NULL;
  int mol_n = 0;
  if (n_aggregates == 0){
    *out_mols = NULL;
    *out_n = 0;
    return 0;
  }
  if (n_aggregates > (size_t)INT_MAX) return -1;
  candidates = (molecule_candidate_t*)malloc(sizeof(*candidates) * n_aggregates);
  mols = (molecule_t*)calloc(n_aggregates, sizeof(*mols));
  if (!candidates || !mols){
    free(candidates);
    free(mols);
    return -1;
  }
  for (size_t i=0; i<n_aggregates; ++i){
    const map_item_t* item = find_in_key_map(km, km_n,
                                             aggregates[i].key,
                                             aggregates[i].umi);
    candidates[i].key = aggregates[i].key;
    if (!item){
      free(candidates);
      free(mols);
      return -1;
    }
    candidates[i].corr = item->corr_umi;
    candidates[i].count = aggregates[i].count;
    candidates[i].best_ordinal = aggregates[i].best_ordinal;
    candidates[i].best_mapq = aggregates[i].best_mapq;
    candidates[i].best_ref_span = aggregates[i].best_ref_span;
    candidates[i].best_query_len = aggregates[i].best_query_len;
  }
  qsort(candidates, n_aggregates, sizeof(*candidates), cmp_molecule_candidate);
  for (size_t i=0; i<n_aggregates; ){
    size_t j = i + 1;
    molecule_t* mol = &mols[mol_n];
    while (j<n_aggregates &&
           strcmp(candidates[j].key, candidates[i].key) == 0 &&
           strcmp(candidates[j].corr, candidates[i].corr) == 0) j++;
    mol->key = sdup(candidates[i].key);
    mol->corr = sdup(candidates[i].corr);
    mol->best_ordinal = -1;
    if (!mol->key || !mol->corr){
      free(candidates);
      free_molecules(mols, mol_n + 1);
      return -1;
    }
    for (size_t k=i; k<j; ++k){
      if (candidates[k].count > INT_MAX - mol->count){
        free(candidates);
        free_molecules(mols, mol_n + 1);
        return -1;
      }
      mol->count += candidates[k].count;
      if (candidate_is_better(&candidates[k], mol)){
        mol->best_ordinal = candidates[k].best_ordinal;
        mol->best_mapq = candidates[k].best_mapq;
        mol->best_ref_span = candidates[k].best_ref_span;
        mol->best_query_len = candidates[k].best_query_len;
      }
    }
    mol_n++;
    i = j;
  }
  free(candidates);
  *out_mols = mols;
  *out_n = mol_n;
  return 0;
}

typedef struct {
  char* read_key;
  char* corr;
  char* molecule_id;
  int is_dup;
} read_status_t;

static void free_read_statuses(read_status_t* statuses, int n){
  if (!statuses) return;
  for (int i=0; i<n; ++i){
    free(statuses[i].read_key);
    free(statuses[i].corr);
    free(statuses[i].molecule_id);
  }
  free(statuses);
}

static int read_pair_role(const bam1_t* b){
  int role = 0;
  if (b->core.flag & BAM_FREAD1) role |= 1;
  if (b->core.flag & BAM_FREAD2) role |= 2;
  return role;
}

static int build_read_identity_buf(const cli_opts_t* o, bam1_t* b,
                                   kstring_t* out){
  const char* qname = bam_get_qname(b);
  const char* scope = "";
  size_t scope_len;
  int role = read_pair_role(b);
  if (!out) return -1;
  if (o->isolate_inputs){
    scope = get_tag_Z(b, o->input_scope_tag);
    if (!scope) return -1;
  }
  if (!qname) qname = "";
  scope_len = strlen(scope);
  out->l = 0;
  return ksprintf(out, "%zu:%s|%d|%s", scope_len, scope, role, qname) < 0 ? -1 : 0;
}

static int cmp_read_status(const void* a, const void* b){
  const read_status_t* x = (const read_status_t*)a;
  const read_status_t* y = (const read_status_t*)b;
  return strcmp(x->read_key, y->read_key);
}

static int read_status_same_assignment(const read_status_t* a, const read_status_t* b){
  if (a->is_dup != b->is_dup || strcmp(a->corr, b->corr) != 0) return 0;
  if (!a->molecule_id && !b->molecule_id) return 1;
  if (!a->molecule_id || !b->molecule_id) return 0;
  return strcmp(a->molecule_id, b->molecule_id) == 0;
}

static int append_read_status(read_status_t** statuses, int* n, int* cap,
                              char* read_key, const char* corr,
                              char* molecule_id, int is_dup){
  read_status_t* grown;
  char* corr_copy;
  if (*n == *cap){
    int new_cap = *cap ? *cap * 2 : 1024;
    grown = (read_status_t*)realloc(*statuses, sizeof(*grown) * (size_t)new_cap);
    if (!grown) return -1;
    *statuses = grown;
    *cap = new_cap;
  }
  corr_copy = sdup(corr);
  if (!corr_copy) return -1;
  (*statuses)[*n].read_key = read_key;
  (*statuses)[*n].corr = corr_copy;
  (*statuses)[*n].molecule_id = molecule_id;
  (*statuses)[*n].is_dup = is_dup;
  (*n)++;
  return 0;
}

static int cmp_string_ptr(const void* a, const void* b){
  const char* const* x = (const char* const*)a;
  const char* const* y = (const char* const*)b;
  return strcmp(*x, *y);
}

static void sort_unique_strvec(strvec_t* values){
  size_t out;
  if (!values || values->n < 2) return;
  qsort(values->data, values->n, sizeof(*values->data), cmp_string_ptr);
  out = 1;
  for (size_t i=1; i<values->n; ++i){
    if (strcmp(values->data[i], values->data[out - 1]) == 0){
      free(values->data[i]);
    } else {
      values->data[out++] = values->data[i];
    }
  }
  values->n = out;
}

static int sorted_strvec_contains(const strvec_t* values, const char* value){
  size_t lo = 0, hi = values->n;
  while (lo < hi){
    size_t mid = lo + (hi - lo) / 2;
    int cmp = strcmp(values->data[mid], value);
    if (cmp == 0) return 1;
    if (cmp < 0) lo = mid + 1;
    else hi = mid;
  }
  return 0;
}

static int build_primary_read_statuses(const cli_opts_t* o, const char* bucket_bam,
                                       const key_map_t* km, int km_n,
                                       molecule_t* mols, int mol_n, int bucket_id,
                                       const strvec_t* non_primary_keys,
                                       read_status_t** out_statuses, int* out_n){
  io_ctx_t io = {0};
  bam1_t* b = NULL;
  read_status_t* statuses = NULL;
  int status_n = 0, status_cap = 0;
  long ordinal = 0;
  int read_rc = 0;
  kstring_t read_identity = {0, 0, NULL};
  if (!non_primary_keys || non_primary_keys->n == 0){
    *out_statuses = NULL;
    *out_n = 0;
    return 0;
  }

  if (io_open(bucket_bam, NULL, &io) != 0) goto fail;
  b = bam_init1();
  if (!b) goto fail;
  read_rc = 0;
  while ((read_rc = sam_read1(io.in, io.hdr, b)) >= 0){
    long current_ordinal = ordinal++;
    const char* raw = get_tag_Z(b, o->umi_tag);
    const char* cb = get_tag_Z(b, o->cell_tag);
    if (!is_unmapped(b) && !is_non_primary(b) && raw && cb){
      char* read_key;
      if (build_read_identity_buf(o, b, &read_identity) != 0) goto fail;
      if (!sorted_strvec_contains(non_primary_keys, read_identity.s)) continue;
      read_key = sdup(read_identity.s);
      if (!read_key) goto fail;
      char* key = build_group_key(b, o->cell_tag, o->gene_tag,
                                  o->source_tag, active_input_scope_tag(o),
                                  o->no_gene, o->no_structure,
                                  o->locus_bin, o->sj_jitter, o->end_bin);
      const map_item_t* item;
      const char* corr;
      int molecule_idx;
      int is_dup;
      char* molecule_id = NULL;
      if (!key){ free(read_key); goto fail; }
      item = find_in_key_map(km, km_n, key, raw);
      if (!item){
        LOG_BUCKET(bucket_id, "primary record %s has no UMI correction mapping",
                   bam_get_qname(b));
        free(read_key);
        free(key);
        goto fail;
      }
      corr = item->corr_umi;
      molecule_idx = molecule_find(mols, mol_n, key, corr);
      if (molecule_idx < 0){
        LOG_BUCKET(bucket_id, "primary record %s has no molecule assignment", bam_get_qname(b));
        free(read_key);
        free(key);
        goto fail;
      }
      is_dup = mols[molecule_idx].best_ordinal != current_ordinal;
      if (o->mol_tag){
        molecule_id = build_molecule_id(cb, key, corr);
        if (!molecule_id){ free(read_key); free(key); goto fail; }
      }
      if (append_read_status(&statuses, &status_n, &status_cap, read_key,
                             corr, molecule_id, is_dup) != 0){
        free(read_key);
        free(molecule_id);
        free(key);
        goto fail;
      }
      free(key);
    }
  }
  if (read_rc < -1) goto fail;
  bam_destroy1(b);
  b = NULL;
  io_close(&io);

  if (status_n > 1){
    qsort(statuses, (size_t)status_n, sizeof(*statuses), cmp_read_status);
    for (int i=0; i<status_n; ){
      int j = i + 1;
      while (j < status_n && strcmp(statuses[j].read_key, statuses[i].read_key) == 0) j++;
      for (int k=i+1; k<j; ++k){
        if (!read_status_same_assignment(&statuses[i], &statuses[k])){
          LOG_BUCKET(bucket_id,
                     "ambiguous primary records for read identity %s; use unique QNAMEs or --isolate-inputs",
                     statuses[i].read_key);
          goto fail;
        }
      }
      i = j;
    }

    int write_idx = 0;
    for (int i=0; i<status_n; ){
      int j = i + 1;
      while (j < status_n && strcmp(statuses[j].read_key, statuses[i].read_key) == 0) j++;
      if (write_idx != i) statuses[write_idx] = statuses[i];
      write_idx++;
      for (int k=i+1; k<j; ++k){
        free(statuses[k].read_key);
        free(statuses[k].corr);
        free(statuses[k].molecule_id);
      }
      i = j;
    }
    status_n = write_idx;
  }
  *out_statuses = statuses;
  *out_n = status_n;
  free(read_identity.s);
  return 0;

fail:
  free(read_identity.s);
  free_read_statuses(statuses, status_n);
  bam_destroy1(b);
  io_close(&io);
  return -1;
}

static const read_status_t* find_read_status(const read_status_t* statuses, int n,
                                             const char* read_key){
  int lo = 0, hi = n - 1;
  while (lo <= hi){
    int mid = lo + (hi - lo) / 2;
    int cmp = strcmp(statuses[mid].read_key, read_key);
    if (cmp == 0) return &statuses[mid];
    if (cmp < 0) lo = mid + 1;
    else hi = mid - 1;
  }
  return NULL;
}

static int phase_bucket_dedup(const cli_opts_t* o){
  int status = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(o->threads) reduction(|:status)
#endif
  for (int i=0;i<o->buckets;i++){
    char ib[4096], ob[4096], molf[4096], asnf[4096], corrf[4096];
    snprintf(ib, sizeof(ib), "%s/bucket_%05d.bam", o->tmp_dir, i);
    snprintf(ob, sizeof(ob), "%s/bucket_%05d.dedup.bam", o->tmp_dir, i);
    snprintf(molf, sizeof(molf), "%s/bucket_%05d.molecules.tsv", o->tmp_dir, i);
    snprintf(asnf, sizeof(asnf), "%s/bucket_%05d.assignments.tsv", o->tmp_dir, i);
    snprintf(corrf, sizeof(corrf), "%s/bucket_%05d.corrections.tsv", o->tmp_dir, i);

    key_map_t* km=NULL; int km_n=0;
    molecule_t* mols=NULL; int mol_n=0;
    aggregate_t* aggregates=NULL; size_t n_aggregates=0;
    read_status_t* read_statuses=NULL; int read_status_n=0;
    strvec_t non_primary_keys;
    int bucket_failed = 0;
    FILE* f_corr = NULL;
    strvec_init(&non_primary_keys);
    if (o->emit_explain){
      f_corr = fopen(corrf, "w");
      if (!f_corr){
        LOG_BUCKET(i, "cannot create corrections TSV");
        strvec_free(&non_primary_keys);
        status |= 1;
        continue;
      }
      fprintf(f_corr, "key\traw_umi\tcorr_umi\traw_count\tseed_count\thamming\traw_avgq\tseed_avgq\tconfidence\treason\tbucket\tmethod\tparent_umi\tparent_count\tedge_hamming\tpath_length\n");
    }
    if (build_bucket_mapping(o, ib, f_corr, i, &km, &km_n,
                             &aggregates, &n_aggregates,
                             &non_primary_keys)!=0){
      if (f_corr){ fclose(f_corr); remove(corrf); }
      strvec_free(&non_primary_keys);
      LOG_BUCKET(i, "mapping failed");
      status |= 1;
      continue;
    }
    if (f_corr){
      int corr_failed = ferror(f_corr);
      if (fclose(f_corr) != 0) corr_failed = 1;
      f_corr = NULL;
      if (corr_failed){
        LOG_BUCKET(i, "failed to flush/close corrections TSV");
        remove(corrf);
        free_aggregates(aggregates, n_aggregates);
        strvec_free(&non_primary_keys);
        free_key_maps(km, km_n);
        status |= 1;
        continue;
      }
    }

    if (build_bucket_molecules_from_aggregates(km, km_n,
                                               aggregates, n_aggregates,
                                               &mols, &mol_n) != 0){
      LOG_BUCKET(i, "molecule representative selection failed");
      if (o->emit_explain){ remove(corrf); }
      free_aggregates(aggregates, n_aggregates);
      strvec_free(&non_primary_keys);
      free_key_maps(km, km_n);
      status |= 1;
      continue;
    }
    free_aggregates(aggregates, n_aggregates);
    aggregates = NULL;
    n_aggregates = 0;

    sort_unique_strvec(&non_primary_keys);
    if (build_primary_read_statuses(o, ib, km, km_n, mols, mol_n, i,
                                    &non_primary_keys,
                                    &read_statuses, &read_status_n) != 0){
      LOG_BUCKET(i, "primary-read status construction failed");
      if (o->emit_explain){ remove(corrf); }
      strvec_free(&non_primary_keys);
      free_molecules(mols, mol_n);
      free_key_maps(km, km_n);
      status |= 1;
      continue;
    }
    strvec_free(&non_primary_keys);

    io_ctx_t io = {0};
    if (io_open_mode(ib, ob, "wb1", &io) != 0){
      LOG_BUCKET(i, "open failed");
      if (o->emit_explain){ remove(corrf); }
      free_read_statuses(read_statuses, read_status_n);
      free_molecules(mols, mol_n);
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
        if (fprintf(f_mol,"key\tumi_corr\tcount\tbucket\n") < 0 ||
            fprintf(f_asn,"qname\tmolecule_id\tdup\tbucket\n") < 0){
          LOG_BUCKET(i, "failed to write TSV headers");
          bucket_failed = 1;
        }
      }
    }

    bam1_t* b = bam_init1();
    if (!b){
      if (f_mol) fclose(f_mol);
      if (f_asn) fclose(f_asn);
      if (o->emit_explain){ remove(corrf); }
      io_close(&io);
      free_read_statuses(read_statuses, read_status_n);
      free_molecules(mols, mol_n);
      free_key_maps(km, km_n);
      status |= 1;
      LOG_BUCKET(i, "out of memory allocating BAM record");
      continue;
    }
    int read_rc = 0;
    long ordinal = 0;
    kstring_t read_identity = {0, 0, NULL};

    while (!bucket_failed && (read_rc = sam_read1(io.in, io.hdr, b)) >= 0){
      long current_ordinal = ordinal++;
      const char* raw = get_tag_Z(b, o->umi_tag);
      const char* cb  = get_tag_Z(b, o->cell_tag);
      const map_item_t* item = NULL;
      const char* corr = NULL;
      if (is_unmapped(b)) {
        if ((raw ? set_tag_Z(b, o->umi_out, raw)
                 : clear_tag_if_present(b, o->umi_out)) < 0 ||
            clear_molecule_tag(o, b) < 0 ||
            set_duplicate_status(o, b, 0) < 0){
          LOG_BUCKET(i, "failed to set unmapped passthrough status");
          bucket_failed = 1;
          break;
        }
        if (clear_input_scope_tag(o, b) < 0){ LOG_BUCKET(i, "failed to remove input scope tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      if (is_non_primary(b)) {
        const read_status_t* inherited = NULL;
        if (build_read_identity_buf(o, b, &read_identity) != 0){ LOG_BUCKET(i, "failed to build non-primary read identity"); bucket_failed = 1; break; }
        inherited = find_read_status(read_statuses, read_status_n, read_identity.s);
        if (inherited){
          if (set_tag_Z(b, o->umi_out, inherited->corr) < 0 ||
              set_duplicate_status(o, b, inherited->is_dup) < 0 ||
              (o->mol_tag && (!inherited->molecule_id ||
                              set_tag_Z(b, o->mol_tag, inherited->molecule_id) < 0))){
            LOG_BUCKET(i, "failed to propagate primary-read status to non-primary record");
            bucket_failed = 1;
            break;
          }
        } else {
          if ((raw ? set_tag_Z(b, o->umi_out, raw)
                   : clear_tag_if_present(b, o->umi_out)) < 0 ||
              clear_molecule_tag(o, b) < 0 ||
              set_duplicate_status(o, b, 0) < 0){
            LOG_BUCKET(i, "failed to set passthrough status on orphan non-primary record");
            bucket_failed = 1;
            break;
          }
        }
        if (clear_input_scope_tag(o, b) < 0){ LOG_BUCKET(i, "failed to remove input scope tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      if (!raw) {
        if (clear_tag_if_present(b, o->umi_out) < 0 ||
            clear_molecule_tag(o, b) < 0 ||
            set_duplicate_status(o, b, 0) < 0){
          LOG_BUCKET(i, "failed to clear derived tags from record without raw UMI");
          bucket_failed = 1;
          break;
        }
        if (clear_input_scope_tag(o, b) < 0){ LOG_BUCKET(i, "failed to remove input scope tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      if (!cb) {
        if (set_tag_Z(b, o->umi_out, raw) < 0 ||
            clear_molecule_tag(o, b) < 0 ||
            set_duplicate_status(o, b, 0) < 0){
          LOG_BUCKET(i, "failed to set passthrough status on record without cell barcode");
          bucket_failed = 1;
          break;
        }
        if (clear_input_scope_tag(o, b) < 0){ LOG_BUCKET(i, "failed to remove input scope tag"); bucket_failed = 1; break; }
        if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); bucket_failed = 1; break; }
        continue;
      }
      char* cb_copy = NULL;
      if (o->mol_tag || o->emit_tsv){
        cb_copy = sdup(cb);
        if (!cb_copy){ LOG_BUCKET(i, "failed to copy cell barcode"); bucket_failed = 1; break; }
      }
      char* key = build_group_key(b, o->cell_tag, o->gene_tag,
                                  o->source_tag, active_input_scope_tag(o),
                                  o->no_gene, o->no_structure,
                                  o->locus_bin, o->sj_jitter, o->end_bin);
      if (!key){ LOG_BUCKET(i, "key build failed"); free(cb_copy); bucket_failed = 1; break; }

      item = find_in_key_map(km, km_n, key, raw);
      if (!item){
        LOG_BUCKET(i, "record %s has no UMI correction mapping", bam_get_qname(b));
        free(cb_copy);
        free(key);
        bucket_failed = 1;
        break;
      }
      corr = item->corr_umi;
      if (corr && set_tag_Z(b, o->umi_out, corr) < 0){ LOG_BUCKET(i, "failed to set corrected UMI tag"); free(cb_copy); free(key); bucket_failed = 1; break; }

      int idx = molecule_find(mols, mol_n, key, corr);
      if (idx < 0){
        LOG_BUCKET(i, "record %s has no molecule assignment", bam_get_qname(b));
        free(cb_copy);
        free(key);
        bucket_failed = 1;
        break;
      }
      int is_dup = mols[idx].best_ordinal != current_ordinal;

      if (set_duplicate_status(o, b, is_dup) < 0){ LOG_BUCKET(i, "failed to set duplicate tag"); free(cb_copy); free(key); bucket_failed = 1; break; }
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
        if (fprintf(f_asn, "%s\t%s\t%d\t%d\n", qn?qn:"*", mi2, is_dup ? 1 : 0, i) < 0){
          free(mi2);
          free(cb_copy);
          free(key);
          LOG_BUCKET(i, "failed to write assignments TSV");
          bucket_failed = 1;
          break;
        }
        free(mi2);
      }

      if (clear_input_scope_tag(o, b) < 0){ LOG_BUCKET(i, "failed to remove input scope tag"); free(cb_copy); free(key); bucket_failed = 1; break; }
      if (sam_write1(io.out, io.hdr, b) < 0){ LOG_BUCKET(i, "write fail"); free(cb_copy); free(key); bucket_failed = 1; break; }
      free(cb_copy);
      free(key);
    }
    if (read_rc < -1){ LOG_BUCKET(i, "read fail"); bucket_failed = 1; }
    free(read_identity.s);
    bam_destroy1(b);
    if (io_close(&io) != 0){
      LOG_BUCKET(i, "failed to flush/close deduplicated BAM");
      bucket_failed = 1;
    }

    if (!bucket_failed && o->emit_tsv && f_mol){
      for (int t=0;t<mol_n;++t){
        if (fprintf(f_mol, "%s\t%s\t%d\t%d\n", mols[t].key, mols[t].corr, mols[t].count, i) < 0){
          LOG_BUCKET(i, "failed to write molecules TSV");
          bucket_failed = 1;
          break;
        }
      }
    }
    if (f_mol && fclose(f_mol) != 0){
      LOG_BUCKET(i, "failed to flush/close molecules TSV");
      bucket_failed = 1;
    }
    if (f_asn && fclose(f_asn) != 0){
      LOG_BUCKET(i, "failed to flush/close assignments TSV");
      bucket_failed = 1;
    }

    free_key_maps(km, km_n);
    free_molecules(mols, mol_n);
    free_read_statuses(read_statuses, read_status_n);

    if (bucket_failed){
      remove(ob);
      if (o->emit_tsv){ remove(molf); remove(asnf); }
      if (o->emit_explain){ remove(corrf); }
      status |= 1;
      LOG_BUCKET(i, "failed");
    } else {
      if (!o->keep_tmp && remove(ib) != 0 && errno != ENOENT){
        LOG_BUCKET(i, "warning: cannot remove completed source bucket %s", ib);
      }
      LOG_BUCKET(i, "done");
    }
  }
  return status ? -1 : 0;
}

static int phase_concat(const cli_opts_t* o){
  char* outp = path_with_suffix(o->out_prefix, ".dedup.bam");
  char* stage = NULL;
  samFile* out = NULL;
  samFile* hin = NULL;
  bam_hdr_t* hdr = NULL;
  bam1_t* b = NULL;
  int failed = 0;
  char hb[4096];
  if (!outp || create_staging_path(outp, &stage) != 0){
    LOG_PIPELINE("[concat] cannot allocate staging output");
    free(outp);
    return -1;
  }
  out = sam_open(stage, "wb");
  if (!out){
    LOG_PIPELINE("[concat] cannot create %s", stage);
    failed = 1;
    goto cleanup;
  }
  if (o->threads > 1 && hts_set_threads(out, o->threads) < 0)
    LOG_PIPELINE("[concat] warning: cannot enable threaded output compression");

  snprintf(hb, sizeof(hb), "%s/bucket_%05d.dedup.bam", o->tmp_dir, 0);
  hin = sam_open(hb, "r");
  if (!hin){
    LOG_PIPELINE("[concat] cannot open header %s", hb);
    failed = 1;
    goto cleanup;
  }
  hdr = sam_hdr_read(hin);
  if (!hdr){
    LOG_PIPELINE("[concat] cannot read header");
    failed = 1;
    goto cleanup;
  }
  if (sam_close(hin) < 0){
    LOG_PIPELINE("[concat] cannot close header source %s", hb);
    failed = 1;
    hin = NULL;
    goto cleanup;
  }
  hin = NULL;
  if (mark_header_sort_unknown(hdr) != 0 || sam_hdr_write(out, hdr) < 0){
    LOG_PIPELINE("[concat] header write failed");
    failed = 1;
    goto cleanup;
  }
  b = bam_init1();
  if (!b){
    failed = 1;
    goto cleanup;
  }
  for (int i=0; i<o->buckets; ++i){
    char ib[4096];
    samFile* in;
    bam_hdr_t* h2;
    int read_rc = 0;
    snprintf(ib, sizeof(ib), "%s/bucket_%05d.dedup.bam", o->tmp_dir, i);
    in = sam_open(ib, "r");
    if (!in){ LOG_PIPELINE("[concat] cannot open %s", ib); failed = 1; break; }
    h2 = sam_hdr_read(in);
    if (!h2 || !headers_compatible(hdr, h2)){
      LOG_PIPELINE("[concat] incompatible or unreadable header in %s", ib);
      bam_hdr_destroy(h2);
      sam_close(in);
      failed = 1;
      break;
    }
    while ((read_rc = sam_read1(in, h2, b)) >= 0){
      if (sam_write1(out, hdr, b) < 0){ LOG_PIPELINE("[concat] write fail"); failed = 1; break; }
    }
    bam_hdr_destroy(h2);
    if (read_rc < -1){ LOG_PIPELINE("[concat] read fail %s", ib); failed = 1; }
    if (sam_close(in) < 0){ LOG_PIPELINE("[concat] close fail %s", ib); failed = 1; }
    if (!failed && !o->keep_tmp && remove(ib) != 0 && errno != ENOENT){
      LOG_PIPELINE("[concat] warning: cannot remove completed deduplicated bucket %s", ib);
    }
    if (failed) break;
  }

cleanup:
  if (hin && sam_close(hin) < 0) failed = 1;
  bam_destroy1(b);
  bam_hdr_destroy(hdr);
  if (out && sam_close(out) < 0){
    LOG_PIPELINE("[concat] failed to flush/close %s", stage);
    failed = 1;
  }
  if (!failed && rename(stage, outp) != 0){
    LOG_PIPELINE("[concat] cannot atomically replace %s: %s", outp, strerror(errno));
    failed = 1;
  }
  if (failed) remove(stage);
  free(stage);
  free(outp);
  return failed ? -1 : 0;
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

static int concat_bucket_report(const cli_opts_t* o, const char* suffix,
                                const char* header_prefix, const char* output_name){
  char* outp = NULL;
  char* stage = NULL;
  FILE* out = NULL;
  FILE* in = NULL;
  int failed = 0;
  int need = snprintf(NULL, 0, "%s.%s", o->out_prefix, output_name);
  if (need < 0) return -1;
  outp = (char*)malloc((size_t)need + 1);
  if (!outp) return -1;
  snprintf(outp, (size_t)need + 1, "%s.%s", o->out_prefix, output_name);
  if (create_staging_path(outp, &stage) != 0){
    free(outp);
    return -1;
  }
  out = fopen(stage, "w");
  if (!out){
    LOG_PIPELINE("[concat] cannot create %s", stage);
    failed = 1;
    goto cleanup;
  }
  for (int i=0; i<o->buckets && !failed; ++i){
    char path[4096];
    char buf[1<<15];
    snprintf(path, sizeof(path), "%s/bucket_%05d.%s", o->tmp_dir, i, suffix);
    in = fopen(path, "r");
    if (!in){
      LOG_PIPELINE("[concat] cannot open %s", path);
      failed = 1;
      break;
    }
    if (fgets(buf,sizeof(buf),in)){
      if ((i == 0 || strncmp(buf, header_prefix, strlen(header_prefix)) != 0) &&
          fputs(buf,out) == EOF){
        LOG_PIPELINE("[concat] write fail %s", outp);
        failed = 1;
      }
    }
    while (!failed && fgets(buf,sizeof(buf),in)){
      if (fputs(buf,out) == EOF){
        LOG_PIPELINE("[concat] write fail %s", outp);
        failed = 1;
        break;
      }
    }
    if (!failed && ferror(in)){
      LOG_PIPELINE("[concat] read fail %s", path);
      failed = 1;
    }
    if (fclose(in) != 0){
      LOG_PIPELINE("[concat] close fail %s", path);
      failed = 1;
    }
    in = NULL;
    if (!failed && !o->keep_tmp && remove(path) != 0 && errno != ENOENT){
      LOG_PIPELINE("[concat] warning: cannot remove completed report bucket %s", path);
    }
  }
  if (out && fclose(out) != 0){
    LOG_PIPELINE("[concat] failed to flush/close %s", outp);
    failed = 1;
  }
  out = NULL;
  if (!failed && rename(stage, outp) != 0){
    LOG_PIPELINE("[concat] cannot atomically replace %s: %s", outp, strerror(errno));
    failed = 1;
  }

cleanup:
  if (in && fclose(in) != 0) failed = 1;
  if (out && fclose(out) != 0) failed = 1;
  if (failed) remove(stage);
  free(stage);
  free(outp);
  return failed ? -1 : 0;
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
  int result = -1;
  if (validate_output_aliases(o) != 0) return -1;
  LOG_PIPELINE("Phase 1: split into %d buckets at %s", o->buckets, o->tmp_dir);
  if (phase_split(o)!=0) return -1;
#ifndef _OPENMP
  if (o->threads > 1){
    LOG_PIPELINE("OpenMP support is not enabled in this build; bucket processing will run serially despite --threads=%d", o->threads);
  }
#endif
  LOG_PIPELINE("Phase 2: per-bucket dedup (threads=%d)", o->threads);
  if (phase_bucket_dedup(o)!=0) goto cleanup;
  LOG_PIPELINE("Phase 3: concat into %s.dedup.bam", o->out_prefix);
  if (phase_concat(o)!=0) goto cleanup;
  if (o->emit_tsv || o->emit_explain){
    LOG_PIPELINE("Phase 4: concat reports for %s", o->out_prefix);
    if (concat_reports(o)!=0) goto cleanup;
  }
  result = 0;

cleanup:
  if (!o->keep_tmp){
    LOG_PIPELINE("Cleanup: removing temporary files under %s%s", o->tmp_dir,
                 result == 0 ? "" : " after failure");
    cleanup_tmp(o);
  }
  if (result == 0) LOG_PIPELINE("DONE");
  return result;
}
