#include "io.h"
#include <stdio.h>

int io_open_mode(const char* inbam, const char* outbam, const char* out_mode,
                 io_ctx_t* io){
  io->in = sam_open(inbam, "r");
  if (!io->in){
    fprintf(stderr, "[io] cannot open %s\n", inbam);
    return -1;
  }

  io->hdr = sam_hdr_read(io->in);
  if (!io->hdr){
    fprintf(stderr, "[io] cannot read header %s\n", inbam);
    sam_close(io->in);
    return -1;
  }

  io->out = NULL;
  if (outbam){
    io->out = sam_open(outbam, out_mode ? out_mode : "wb");
    if (!io->out){
      fprintf(stderr, "[io] cannot create %s\n", outbam);
      io_close(io);
      return -1;
    }
    if (sam_hdr_write(io->out, io->hdr) < 0){
      fprintf(stderr, "[io] hdr write fail %s\n", outbam);
      io_close(io);
      return -1;
    }
  }
  return 0;
}

int io_open(const char* inbam, const char* outbam, io_ctx_t* io){
  return io_open_mode(inbam, outbam, "wb", io);
}

int io_close(io_ctx_t* io){
  int failed = 0;
  if (!io) return 0;
  if (io->in && sam_close(io->in) < 0) failed = 1;
  if (io->out && sam_close(io->out) < 0) failed = 1;
  if (io->hdr) bam_hdr_destroy(io->hdr);
  io->in = io->out = NULL;
  io->hdr = NULL;
  return failed ? -1 : 0;
}
