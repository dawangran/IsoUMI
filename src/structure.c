#include "structure.h"
#include "key_build.h"
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char* scope;
  const char* neighborhood;
  const char* raw;
  hts_pos_t start, end;
  hts_pos_t* junctions;
  int nj;
  long ordinal;
} observation_t;
typedef struct { size_t first, end; int group; } node_t;
typedef struct {
  int* nodes;
  int n, cap, representative;
  size_t support;
  int* witnesses;
  int nw, cw;
  int* assigned;
  int na, ca;
} anchor_t;

enum { ATTACHED = 0, AMBIGUOUS = 1, UNSUPPORTED = 2 };

static int cmp_pos(hts_pos_t a, hts_pos_t b){ return (a > b) - (a < b); }
static int cmp_observation(const void* va, const void* vb){
  const observation_t* a = va;
  const observation_t* b = vb;
  int c = strcmp(a->scope, b->scope);
  if (c) return c;
  c = strcmp(a->neighborhood, b->neighborhood);
  if (c) return c;
  if (a->nj != b->nj) return (b->nj > a->nj) - (b->nj < a->nj);
  c = cmp_pos(b->end - b->start, a->end - a->start);
  if (c) return c;
  c = cmp_pos(a->start, b->start);
  if (c) return c;
  c = cmp_pos(a->end, b->end);
  if (c) return c;
  for (int i=0; i<2*a->nj; ++i){
    c = cmp_pos(a->junctions[i], b->junctions[i]);
    if (c) return c;
  }
  c = strcmp(a->raw, b->raw);
  if (c) return c;
  return (a->ordinal > b->ordinal) - (a->ordinal < b->ordinal);
}
static int same_geometry(const observation_t* a, const observation_t* b){
  return a->start == b->start && a->end == b->end && a->nj == b->nj &&
    (!a->nj || memcmp(a->junctions, b->junctions,
                     (size_t)a->nj * 2 * sizeof(hts_pos_t)) == 0);
}
static int within(hts_pos_t a, hts_pos_t b, int tolerance){
  return a >= b ? a - b <= tolerance : b - a <= tolerance;
}
static int overlap(hts_pos_t a, hts_pos_t b, hts_pos_t c, hts_pos_t d){
  return a < d && c < b;
}
static hts_pos_t exon_start(const observation_t* a, int i){
  return i ? a->junctions[2*i-1] : a->start;
}
static hts_pos_t exon_end(const observation_t* a, int i){
  return i < a->nj ? a->junctions[2*i] : a->end;
}
/* Exonic sequence cannot run into the interior of another read's intron.
   This also rejects a partial retained intron, not just wholly covered introns. */
static int exon_in_intron(const observation_t* a, const observation_t* b, int tol){
  for (int j=0; j<b->nj; ++j){
    hts_pos_t left = b->junctions[2*j], right = b->junctions[2*j+1];
    if (right - left <= (hts_pos_t)2 * tol) continue;
    left += tol; right -= tol;
    for (int i=0; i<=a->nj; ++i)
      if (overlap(exon_start(a,i), exon_end(a,i), left, right)) return 1;
  }
  return 0;
}
static int junction_match(const observation_t* a, int i,
                          const observation_t* b, int j, int tol){
  return within(a->junctions[2*i],b->junctions[2*j],tol) &&
         within(a->junctions[2*i+1],b->junctions[2*j+1],tol);
}
static int junction_observed(const observation_t* a, int i, const observation_t* b){
  return a->junctions[2*i] >= b->start && a->junctions[2*i+1] <= b->end;
}
static int chain_conflict(const observation_t* a, const observation_t* b, int tol){
  int required=0;
  for (int i=0; i<a->nj; ++i) if (junction_observed(a,i,b)) { required=1; break; }
  for (int j=0; j<b->nj && !required; ++j) if (junction_observed(b,j,a)) required=1;
  /* Terminal overlap may contain no complete intron from either read. */
  if (!required) return 0;
  if (!a->nj || !b->nj) return 1;
  /* Try every continuous chain offset, rather than greedily matching a junction
     to the first boundary within tolerance. This remains symmetric and preserves
     one-to-one matching even for introns shorter than the tolerance. */
  for (int offset=1-a->nj; offset<b->nj; ++offset){
    int valid=1, first=-1, last=-1;
    for (int i=0; i<a->nj; ++i){
      int j=i+offset;
      int matched=j>=0 && j<b->nj && junction_match(a,i,b,j,tol);
      if (matched){
        if (last>=0 && i!=last+1) { valid=0; break; }
        if (first<0) first=i;
        last=i;
      } else if (junction_observed(a,i,b)){ valid=0; break; }
    }
    if (!valid || first<0) continue;
    for (int j=0; j<b->nj; ++j){
      int i=j-offset;
      if (junction_observed(b,j,a) &&
          (i<0 || i>=a->nj || !junction_match(a,i,b,j,tol))){ valid=0; break; }
    }
    if (valid) return 0;
  }
  return 1;
}
static int contradicts(const observation_t* a, const observation_t* b, int tol){
  return exon_in_intron(a,b,tol) || exon_in_intron(b,a,tol) || chain_conflict(a,b,tol);
}

static int compatible(const observation_t* a, const observation_t* b, int tol){
  int exonic_overlap=0;
  for (int i=0; i<=a->nj && !exonic_overlap; ++i)
    for (int j=0; j<=b->nj; ++j)
      if (overlap(exon_start(a,i),exon_end(a,i),exon_start(b,j),exon_end(b,j))){
        exonic_overlap=1; break;
      }
  return exonic_overlap && !contradicts(a,b,tol);
}

int structure_map_add(structure_map_t* map, char* scope,
                      const char* neighborhood, const char* raw,
                      const bam1_t* b, long ordinal){
  observation_t value = {0};
  observation_t* obs;
  hts_pos_t pos;
  uint32_t* cigar;
  if (!scope || !neighborhood || !raw || !b || ordinal < 0){ free(scope); return -1; }
  value.scope=scope; value.neighborhood=neighborhood; value.raw=raw;
  value.start=b->core.pos; value.end=bam_endpos(b); value.ordinal=ordinal;
  cigar=bam_get_cigar(b);
  for (uint32_t i=0; i<b->core.n_cigar; ++i)
    if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP){
      if (value.nj == INT_MAX/2){ free(scope); return -1; }
      value.nj++;
    }
  if (value.nj){
    value.junctions=malloc((size_t)value.nj * 2 * sizeof(hts_pos_t));
    if (!value.junctions){ free(scope); return -1; }
  }
  pos=value.start;
  int n=0;
  for (uint32_t i=0; i<b->core.n_cigar; ++i){
    int op=bam_cigar_op(cigar[i]);
    hts_pos_t length=bam_cigar_oplen(cigar[i]);
    if (op == BAM_CREF_SKIP){ value.junctions[n++]=pos; value.junctions[n++]=pos+length; }
    if (bam_cigar_type(op) & 2) pos += length;
  }
  if (map->n == map->cap){
    size_t cap=map->cap ? map->cap*2 : 256;
    if (cap < map->cap || cap > SIZE_MAX/sizeof(*obs)){ free(scope); free(value.junctions); return -1; }
    obs=realloc(map->observations,cap*sizeof(*obs));
    if (!obs){ free(scope); free(value.junctions); return -1; }
    map->observations=obs; map->cap=cap;
  }
  obs=map->observations;
  obs[map->n++]=value;
  if ((size_t)ordinal >= map->ordinal_count) map->ordinal_count=(size_t)ordinal+1;
  return 0;
}
static int append_node(int** nodes, int* n, int* cap, int node){
  if (*n == *cap){
    if (*cap > INT_MAX/2) return -1;
    int next=*cap ? *cap*2 : 8;
    int* grown=realloc(*nodes,(size_t)next*sizeof(int));
    if (!grown) return -1;
    *nodes=grown; *cap=next;
  }
  (*nodes)[(*n)++]=node;
  return 0;
}
static int anchor_add(anchor_t* anchor, int node, size_t support){
  if (append_node(&anchor->nodes,&anchor->n,&anchor->cap,node)!=0) return -1;
  anchor->support+=support;
  return 0;
}
static int no_conflict(const int* members, int n, const node_t* nodes,
                       const observation_t* obs, const observation_t* candidate, int tol){
  for (int j=0; j<n; ++j)
    if (contradicts(&obs[nodes[members[j]].first],candidate,tol)) return 0;
  return 1;
}
static int anchor_compatible(const anchor_t* anchor, const node_t* nodes,
                            const observation_t* obs, const observation_t* candidate, int tol){
  /* Fixed geometry supplies positive overlap; witnesses veto only explicit
     conflicts. Opposite-end fragments without overlap are not contradictory.
     Discovery witnesses are frozen, even when subsequently ambiguous, and
     are distinct from the final uniquely assigned support used for ranking. */
  return compatible(&obs[nodes[anchor->representative].first],candidate,tol) &&
    no_conflict(anchor->witnesses,anchor->nw,nodes,obs,candidate,tol);
}
static char* group_key(const observation_t* first, int group){
  char* escaped=escape_key_component(first->neighborhood);
  char* key;
  int need;
  if (!escaped) return NULL;
  need=snprintf(NULL,0,"%s|CP=%s|SG=%d",first->scope,escaped,group);
  if (need < 0){ free(escaped); return NULL; }
  key=malloc((size_t)need+1);
  if (key) snprintf(key,(size_t)need+1,"%s|CP=%s|SG=%d",first->scope,escaped,group);
  free(escaped); return key;
}
static size_t raw_support(const anchor_t* anchor, const node_t* nodes,
                          const observation_t* obs, const char* raw){
  size_t count=0;
  for (int i=0; i<anchor->n; ++i){
    const node_t* node=&nodes[anchor->nodes[i]];
    for (size_t j=node->first; j<node->end; ++j)
      if (strcmp(obs[j].raw,raw)==0) count++;
  }
  return count;
}
static int family_finish(structure_map_t* map, size_t first, size_t end, int tol, int min_support){
  observation_t* obs=map->observations;
  if (end-first > (size_t)INT_MAX || end-first > SIZE_MAX/sizeof(anchor_t) ||
      end-first > SIZE_MAX/sizeof(node_t)) return -1;
  node_t* nodes=calloc(end-first,sizeof(*nodes));
  anchor_t* anchors=calloc(end-first,sizeof(*anchors));
  char** groupkeys=NULL;
  int nn=0, na=0, result=-1;
  if (!nodes || !anchors) goto done;
  for (size_t i=first; i<end; ){
    size_t j=i+1;
    while (j<end && same_geometry(&obs[i],&obs[j])) j++;
    nodes[nn++]=(node_t){i,j,-1}; i=j;
  }
  for (int i=0; i<nn; ++i){
    int chosen=-1, matches=0;
    for (int a=0; a<na; ++a)
      if (anchor_compatible(&anchors[a],nodes,obs,&obs[nodes[i].first],tol)) { chosen=a; matches++; }
    if (!matches) { chosen=na++; anchors[chosen].representative=i; }
    if (matches <= 1 && append_node(&anchors[chosen].witnesses,
          &anchors[chosen].nw,&anchors[chosen].cw,i)!=0) goto done;
  }
  /* A later anchor can make an earlier unique attachment ambiguous. Recheck
     every node against the complete fixed anchor set before counting support. */
  for (int i=0; i<nn; ++i){
    int chosen=-1, matches=0;
    for (int a=0; a<na; ++a)
      if (anchor_compatible(&anchors[a],nodes,obs,&obs[nodes[i].first],tol)) { chosen=a; matches++; }
    if (matches != 1) continue;
    /* A node deferred during discovery need not be a witness. Even finally
       unique candidates must respect the evidence already accepted here. */
    if (!no_conflict(anchors[chosen].assigned,anchors[chosen].na,
                     nodes,obs,&obs[nodes[i].first],tol)) continue;
    if (anchor_add(&anchors[chosen],i,nodes[i].end-nodes[i].first)!=0 ||
        append_node(&anchors[chosen].assigned,&anchors[chosen].na,
                    &anchors[chosen].ca,i)!=0) goto done;
    nodes[i].group=chosen;
  }
  groupkeys=calloc((size_t)nn,sizeof(*groupkeys));
  if (!groupkeys) goto done;
  for (int a=0; a<na; ++a){
    groupkeys[a]=group_key(&obs[first],a);
    if (!groupkeys[a]) goto done;
  }
  const int frozen_na=na;
  /* Fix all unique assignments first. Ambiguous decisions never change support,
     but every assignment adds a negative conflict guard for subsequent reads. */
  for (int i=0; i<nn; ++i){
    int chosen=nodes[i].group;
    if (chosen < 0) continue;
    for (size_t j=nodes[i].first; j<nodes[i].end; ++j){
      map->keys[obs[j].ordinal]=groupkeys[chosen];
      map->statuses[obs[j].ordinal]=anchors[chosen].support < (size_t)min_support ? UNSUPPORTED : ATTACHED;
    }
  }
  for (int i=0; i<nn; ++i){
    int fallback=-1;
    if (nodes[i].group >= 0) continue;
    for (size_t j=nodes[i].first; j<nodes[i].end; ++j){
      int chosen=-1;
      size_t best_raw=0, best_total=0;
      int best_qualified=-1;
      unsigned char status=AMBIGUOUS;
      for (int a=0; a<frozen_na; ++a){
        size_t raw;
        int qualified;
        if (!anchor_compatible(&anchors[a],nodes,obs,&obs[j],tol) ||
            !no_conflict(anchors[a].assigned,anchors[a].na,nodes,obs,&obs[j],tol)) continue;
        qualified=anchors[a].support >= (size_t)min_support;
        raw=raw_support(&anchors[a],nodes,obs,obs[j].raw);
        if (chosen < 0 || qualified > best_qualified ||
            (qualified == best_qualified && (raw > best_raw ||
              (raw == best_raw && anchors[a].support > best_total)))){
          chosen=a; best_qualified=qualified; best_raw=raw; best_total=anchors[a].support;
        }
      }
      if (chosen < 0){
        if (fallback < 0){
          if (na >= nn) goto done;
          fallback=na++;
          groupkeys[fallback]=group_key(&obs[first],fallback);
          if (!groupkeys[fallback]) goto done;
        }
        chosen=fallback; status=UNSUPPORTED;
      }
      if (!anchors[chosen].na || anchors[chosen].assigned[anchors[chosen].na-1]!=i)
        if (append_node(&anchors[chosen].assigned,&anchors[chosen].na,&anchors[chosen].ca,i)!=0) goto done;
      map->keys[obs[j].ordinal]=groupkeys[chosen];
      map->statuses[obs[j].ordinal]=status;
    }
  }
  {
    char** grown=realloc(map->owned_keys,(map->group_count+(size_t)na)*sizeof(char*));
    if (!grown) goto done;
    map->owned_keys=grown;
    memcpy(grown+map->group_count,groupkeys,(size_t)na*sizeof(char*));
    map->group_count+=(size_t)na;
    free(groupkeys); groupkeys=NULL;
  }
  result=0;
done:
  if (groupkeys){ for (int a=0; a<na; ++a) free(groupkeys[a]); free(groupkeys); }
  if (anchors){ for (int a=0; a<na; ++a){
    free(anchors[a].nodes); free(anchors[a].witnesses); free(anchors[a].assigned);
  }}
  free(anchors); free(nodes);
  return result;
}
int structure_map_finish(structure_map_t* map, int tolerance, int min_support){
  observation_t* obs=map->observations;
  if (!map->n) return 0;
  map->keys=calloc(map->ordinal_count,sizeof(char*));
  map->statuses=calloc(map->ordinal_count,sizeof(unsigned char));
  if (!map->keys || !map->statuses) return -1;
  qsort(obs,map->n,sizeof(*obs),cmp_observation);
  for (size_t i=0; i<map->n; ){
    size_t j=i+1;
    while (j<map->n && strcmp(obs[i].scope,obs[j].scope)==0 &&
           strcmp(obs[i].neighborhood,obs[j].neighborhood)==0) j++;
    if (family_finish(map,i,j,tolerance,min_support)!=0) return -1;
    i=j;
  }
  for (size_t i=0; i<map->n; ++i){ free(obs[i].scope); free(obs[i].junctions); }
  free(obs); map->observations=NULL; map->n=map->cap=0;
  return 0;
}
const char* structure_map_key(const structure_map_t* map, long ordinal){
  return ordinal>=0 && (size_t)ordinal<map->ordinal_count ? map->keys[ordinal] : NULL;
}
const char* structure_map_status(const structure_map_t* map, long ordinal){
  if (ordinal<0 || (size_t)ordinal>=map->ordinal_count || !map->keys[ordinal]) return "unassigned";
  switch (map->statuses[ordinal]){
    case AMBIGUOUS: return "ambiguous";
    case UNSUPPORTED: return "unsupported";
    default: return "compatible";
  }
}
void structure_map_destroy(structure_map_t* map){
  observation_t* obs=map->observations;
  for (size_t i=0; i<map->n; ++i){ free(obs[i].scope); free(obs[i].junctions); }
  for (size_t i=0; i<map->group_count; ++i) free(map->owned_keys[i]);
  free(obs); free(map->owned_keys); free(map->keys); free(map->statuses);
  memset(map,0,sizeof(*map));
}
