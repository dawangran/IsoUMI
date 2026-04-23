#include "key.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

const char* get_tag_Z(bam1_t* b, const char* tag){
  if (!tag || strlen(tag)!=2) return NULL;
  uint8_t* p = bam_aux_get(b, tag);
  if (!p) return NULL;
  return bam_aux2Z(p);
}
int set_tag_Z(bam1_t* b, const char* tag, const char* val){
  if (!tag || strlen(tag)!=2) return -1;
  if (!val) val="";
  return bam_aux_update_str(b, tag, (int)strlen(val)+1, val);
}
int set_tag_i(bam1_t* b, const char* tag, int32_t val){
  if (!tag || strlen(tag)!=2) return -1;
  return bam_aux_update_int(b, tag, (int64_t)val);
}
uint64_t hash64(const void* data, size_t len){
  const unsigned char* p = (const unsigned char*)data;
  uint64_t h = 1469598103934665603ull;
  for (size_t i=0;i<len;++i){ h ^= p[i]; h *= 1099511628211ull; }
  return h;
}
