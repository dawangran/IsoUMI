#pragma once
#include <htslib/sam.h>
typedef struct { samFile *in, *out; bam_hdr_t *hdr; } io_ctx_t;
int io_open(const char* inbam, const char* outbam, io_ctx_t* io);
void io_close(io_ctx_t* io);
