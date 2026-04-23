#include "vector.h"
#include <stdlib.h>
#include <string.h>

static char* sdup(const char* s){
  size_t n=s?strlen(s):0; char* p=(char*)malloc(n+1);
  if(!p) return NULL; memcpy(p,s?s:"",n+1); return p;
}
void strvec_init(strvec_t *v){ v->data=NULL; v->n=0; v->cap=0; }
int strvec_push(strvec_t *v, const char *s){
  char* dup = sdup(s);
  if(!dup) return -1;
  if (v->n==v->cap){
    size_t nc=v->cap? v->cap*2:8;
    char** nd=(char**)realloc(v->data,sizeof(char*)*nc);
    if(!nd){ free(dup); return -1; } v->data=nd; v->cap=nc;
  }
  v->data[v->n++]=dup; return 0;
}
void strvec_free(strvec_t *v){
  for (size_t i=0;i<v->n;++i) free(v->data[i]);
  free(v->data); v->data=NULL; v->n=v->cap=0;
}
