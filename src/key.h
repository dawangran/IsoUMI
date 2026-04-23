#pragma once
#include <stdint.h>
#include <htslib/sam.h>
const char* get_tag_Z(bam1_t* b, const char* tag);
int  set_tag_Z(bam1_t* b, const char* tag, const char* val);
int  set_tag_i(bam1_t* b, const char* tag, int32_t val);
uint64_t hash64(const void* data, size_t len);
