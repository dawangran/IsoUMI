#pragma once
#include <stdlib.h>
typedef struct { int *p, *sz, n; } uf_t;
static inline int uf_init(uf_t* uf, int n){
  uf->n = n;
  uf->p = NULL;
  uf->sz = NULL;
  if (n <= 0) return 0;
  uf->p = (int*)malloc(sizeof(int) * (size_t)n);
  uf->sz = (int*)malloc(sizeof(int) * (size_t)n);
  if (!uf->p || !uf->sz){
    free(uf->p);
    free(uf->sz);
    uf->p = uf->sz = NULL;
    uf->n = 0;
    return -1;
  }
  for(int i=0;i<n;++i){ uf->p[i]=i; uf->sz[i]=1; }
  return 0;
}
static inline int uf_find(uf_t* uf, int x){
  while(x!=uf->p[x]){ uf->p[x]=uf->p[uf->p[x]]; x=uf->p[x]; }
  return x;
}
static inline void uf_union(uf_t* uf, int a, int b){
  a=uf_find(uf,a); b=uf_find(uf,b); if(a==b) return;
  if(uf->sz[a]<uf->sz[b]){ int t=a;a=b;b=t; }
  uf->p[b]=a; uf->sz[a]+=uf->sz[b];
}
static inline void uf_free(uf_t* uf){
  free(uf->p); free(uf->sz); uf->p=uf->sz=NULL; uf->n=0;
}
