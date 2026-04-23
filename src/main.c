#include "cli.h"
#include "pipeline.h"
int main(int argc, char** argv){
  cli_opts_t o;
  if (parse_args(argc, argv, &o)!=0) return 2;
  int rc = run_pipeline(&o);
  free_opts(&o);
  return rc==0?0:1;
}
