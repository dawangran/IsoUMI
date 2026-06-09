#include "cli.h"
#include "vector.h"
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

static const char* ISOUMI_VERSION = "IsoUMI 0.1.0";

static char* xstrdup(const char* s){
  if(!s) return NULL;
  size_t n=strlen(s);
  char* p=(char*)malloc(n+1);
  if(!p) return NULL;
  memcpy(p,s,n+1);
  return p;
}

static int set_string_opt(char **dst, const char *src){
  char *dup = xstrdup(src);
  if (!dup){
    fprintf(stderr, "out of memory while parsing options\n");
    return -1;
  }
  free(*dst);
  *dst = dup;
  return 0;
}

static int set_tag_opt(char **dst, const char *name, const char *tag){
  if (!tag || strlen(tag) != 2){
    fprintf(stderr, "%s must be a 2-character SAM tag\n", name);
    return -1;
  }
  return set_string_opt(dst, tag);
}

static int parse_int_opt(const char *name, const char *value, int min_value, int *out){
  char *end = NULL;
  long parsed;
  errno = 0;
  parsed = strtol(value, &end, 10);
  if (errno != 0 || end == value || *end != '\0' || parsed < (long)min_value || parsed > INT_MAX){
    fprintf(stderr, "invalid value for %s: %s\n", name, value);
    return -1;
  }
  *out = (int)parsed;
  return 0;
}

static int parse_double_opt(const char *name, const char *value, double min_value, double *out){
  char *end = NULL;
  double parsed;
  errno = 0;
  parsed = strtod(value, &end);
  if (errno != 0 || end == value || *end != '\0' || !isfinite(parsed) || parsed < min_value){
    fprintf(stderr, "invalid value for %s: %s\n", name, value);
    return -1;
  }
  *out = parsed;
  return 0;
}

static int push_bam_path(cli_opts_t *o, const char *path){
  if (strvec_push(&o->bam_list, path) != 0){
    fprintf(stderr, "out of memory while recording BAM input\n");
    return -1;
  }
  return 0;
}

static char* trim_in_place(char *s){
  char *end;
  while (*s && isspace((unsigned char)*s)) s++;
  end = s + strlen(s);
  while (end > s && isspace((unsigned char)end[-1])) --end;
  *end = '\0';
  return s;
}

static void usage(void){
  fprintf(stderr,
"IsoUMI: Isoform-aware long-read UMI correction with cell-bucket sharding\n"
"%s\n"
"\n"
"USAGE:\n"
"  isoumi --bam <in.bam> [--bam <in2.bam> ...] --out <prefix> [OPTIONS]\n"
"  isoumi --bam-list <list.txt> --out <prefix> [OPTIONS]\n"
"\n"
"INPUTS (choose one mode):\n"
"  --bam <FILE>           Input BAM. Repeatable. Example: --bam A.bam --bam B.bam\n"
"  --bam-list <FILE>      Text file; one BAM path per line (supports # comments).\n"
"\n"
"REQUIRED:\n"
"  -o, --out <PREFIX>     Output prefix. Final BAM: <PREFIX>.dedup.bam\n"
"\n"
"PERFORMANCE & SHARDING:\n"
"  --threads <INT>        Worker threads for per-bucket step (default: 4)\n"
"  --buckets <INT>        Number of CB buckets (default: 64)\n"
"  --tmp-dir <DIR>        Bucket directory (default: tmp_buckets)\n"
"\n"
"TAG NAMES (10x-style defaults):\n"
"  --cell-tag <TAG>       Cell barcode tag (default: CB)\n"
"  --umi-tag  <TAG>       Raw UMI tag (default: UR)\n"
"  --umi-qual-tag <TAG>   UMI quality tag (default: UY)\n"
"  --gene-tag <TAG>       Gene tag (default: GX). Use --no-gene to ignore gene\n"
"  --umi-out <TAG>        Corrected UMI tag to write (default: UB)\n"
"  --dup-flag <TAG>       Duplicate flag tag (i32 0/1) (default: DA)\n"
"  --mol-tag <TAG>        Optional; write molecule id (CB|key|UMIcorr)\n"
"\n"
"CORRECTION SETTINGS:\n"
"  --ham <INT>            Max Hamming distance (default: 1)\n"
"  --ratio <FLOAT>        Collapse if smaller/larger <= ratio (default: 0.10)\n"
"  --min-merge-confidence <FLOAT>\n"
"                         Optional confidence floor for quality-aware merges (default: 0.00)\n"
"  --no-quality-aware     Ignore UMI quality when ranking/explaining merges\n"
"  --no-gene              Ignore gene in grouping\n"
"  --locus-bin <INT>      Bin (bp) for non-spliced reads (default: 1000)\n"
"  --sj-jitter <INT>      Round SJ boundaries by this jitter (bp) before hashing (default: 10)\n"
"  --end-bin <INT>        Add strand-aware transcript end bins to grouping key (default: off)\n"
"\n"
"OUTPUTS:\n"
"  --emit-tsv             Also write <out>.molecules.tsv and <out>.assignments.tsv\n"
"  --emit-explain         Also write <out>.corrections.tsv with merge confidence/details\n"
"  --keep-tmp             Keep bucket temp files (default: delete after finish)\n"
"\n"
"HELP:\n"
"  -h, --help             Show this help and exit\n"
"      --version          Print version and exit\n", ISOUMI_VERSION);
}

int parse_args(int argc, char **argv, cli_opts_t *o){
  memset(o, 0, sizeof(*o));
  strvec_init(&o->bam_list);
  o->out_prefix = NULL;
  o->tmp_dir = xstrdup("tmp_buckets");
  o->buckets = 64;
  o->threads = 4;
  o->cell_tag = xstrdup("CB");
  o->umi_tag  = xstrdup("UR");
  o->umi_qual_tag = xstrdup("UY");
  o->gene_tag = xstrdup("GX");
  o->umi_out  = xstrdup("UB");
  o->dup_flag = xstrdup("DA");
  o->mol_tag  = NULL;
  o->no_gene = 0;
  o->ham = 1;
  o->ratio = 0.1;
  o->min_merge_confidence = 0.0;
  o->locus_bin = 1000;
  o->sj_jitter = 10;
  o->end_bin = 0;
  o->emit_tsv = 0;
  o->emit_explain = 0;
  o->quality_aware = 1;
  o->keep_tmp = 0;

  if (!o->tmp_dir || !o->cell_tag || !o->umi_tag || !o->umi_qual_tag || !o->gene_tag || !o->umi_out || !o->dup_flag){
    fprintf(stderr, "out of memory while initializing defaults\n");
    return -1;
  }

  static struct option long_opts[] = {
    {"bam", required_argument, 0, 'b'},
    {"bam-list", required_argument, 0, 1},
    {"out", required_argument, 0, 'o'},
    {"tmp-dir", required_argument, 0, 2},
    {"buckets", required_argument, 0, 3},
    {"threads", required_argument, 0, 4},
    {"cell-tag", required_argument, 0, 5},
    {"umi-tag", required_argument, 0, 6},
    {"umi-qual-tag", required_argument, 0, 7},
    {"gene-tag", required_argument, 0, 8},
    {"umi-out", required_argument, 0, 9},
    {"dup-flag", required_argument, 0, 10},
    {"mol-tag", required_argument, 0, 11},
    {"no-gene", no_argument, 0, 12},
    {"ham", required_argument, 0, 13},
    {"ratio", required_argument, 0, 14},
    {"locus-bin", required_argument, 0, 15},
    {"version", no_argument, 0, 16},
    {"help", no_argument, 0, 'h'},
    {"sj-jitter", required_argument, 0, 17},
    {"emit-tsv", no_argument, 0, 18},
    {"keep-tmp", no_argument, 0, 19},
    {"min-merge-confidence", required_argument, 0, 20},
    {"emit-explain", no_argument, 0, 21},
    {"end-bin", required_argument, 0, 22},
    {"no-quality-aware", no_argument, 0, 23},
    {0,0,0,0}
  };

  int idx=0, c;
  while ((c = getopt_long(argc, argv, "b:o:h", long_opts, &idx)) != -1){
    switch (c){
      case 'b':
        if (push_bam_path(o, optarg) != 0) return -1;
        break;
      case 'o':
        if (set_string_opt(&o->out_prefix, optarg) != 0) return -1;
        break;
      case 'h': usage(); exit(0);
      case 1: { FILE* fp = fopen(optarg, "r");
        if (!fp){ fprintf(stderr,"cannot open bam-list %s\n", optarg); return -1; }
        char buf[8192];
        while (fgets(buf,sizeof(buf),fp)){
          char *nl = strchr(buf,'\n'); if (nl) *nl=0;
          char *line = trim_in_place(buf);
          if (!line[0] || line[0]=='#') continue;
          if (push_bam_path(o, line) != 0){ fclose(fp); return -1; }
        }
        fclose(fp); break; }
      case 2:
        if (set_string_opt(&o->tmp_dir, optarg) != 0) return -1;
        break;
      case 3:
        if (parse_int_opt("--buckets", optarg, 1, &o->buckets) != 0) return -1;
        break;
      case 4:
        if (parse_int_opt("--threads", optarg, 1, &o->threads) != 0) return -1;
        break;
      case 5:
        if (set_tag_opt(&o->cell_tag, "--cell-tag", optarg) != 0) return -1;
        break;
      case 6:
        if (set_tag_opt(&o->umi_tag, "--umi-tag", optarg) != 0) return -1;
        break;
      case 7:
        if (set_tag_opt(&o->umi_qual_tag, "--umi-qual-tag", optarg) != 0) return -1;
        break;
      case 8:
        if (set_tag_opt(&o->gene_tag, "--gene-tag", optarg) != 0) return -1;
        break;
      case 9:
        if (set_tag_opt(&o->umi_out, "--umi-out", optarg) != 0) return -1;
        break;
      case 10:
        if (set_tag_opt(&o->dup_flag, "--dup-flag", optarg) != 0) return -1;
        break;
      case 11:
        if (set_tag_opt(&o->mol_tag, "--mol-tag", optarg) != 0) return -1;
        break;
      case 12:
        o->no_gene=1;
        break;
      case 13:
        if (parse_int_opt("--ham", optarg, 0, &o->ham) != 0) return -1;
        break;
      case 14:
        if (parse_double_opt("--ratio", optarg, 0.0, &o->ratio) != 0) return -1;
        break;
      case 15:
        if (parse_int_opt("--locus-bin", optarg, 1, &o->locus_bin) != 0) return -1;
        break;
      case 16:
        printf("%s\n", ISOUMI_VERSION); exit(0);
      case 17:
        if (parse_int_opt("--sj-jitter", optarg, 0, &o->sj_jitter) != 0) return -1;
        break;
      case 18: o->emit_tsv=1; break;
      case 19: o->keep_tmp=1; break;
      case 20:
        if (parse_double_opt("--min-merge-confidence", optarg, 0.0, &o->min_merge_confidence) != 0) return -1;
        if (o->min_merge_confidence > 1.0){
          fprintf(stderr, "--min-merge-confidence must be <= 1.0\n");
          return -1;
        }
        break;
      case 21: o->emit_explain=1; break;
      case 22:
        if (parse_int_opt("--end-bin", optarg, 1, &o->end_bin) != 0) return -1;
        break;
      case 23: o->quality_aware=0; break;
      default: usage(); return -1;
    }
  }

  if (!o->out_prefix){ fprintf(stderr,"--out is required\n"); usage(); return -1; }
  if (o->bam_list.n==0){ fprintf(stderr,"no BAM inputs (use --bam or --bam-list)\n"); usage(); return -1; }
  return 0;
}

void free_opts(cli_opts_t *o){
  if (!o) return;
  free(o->out_prefix);
  free(o->tmp_dir);
  free(o->cell_tag); free(o->umi_tag); free(o->umi_qual_tag); free(o->gene_tag);
  free(o->umi_out);  free(o->dup_flag);
  free(o->mol_tag);
  strvec_free(&o->bam_list);
}
