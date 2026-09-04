#pragma once
#include <htslib/sam.h>
typedef struct { samFile *in, *out; bam_hdr_t *hdr; } io_ctx_t;
int io_open(const char* inbam, const char* outbam, io_ctx_t* io);
int io_open_mode(const char* inbam, const char* outbam, const char* out_mode,
                 io_ctx_t* io);
int io_close(io_ctx_t* io);
