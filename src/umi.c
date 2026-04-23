#include "umi.h"
#include <string.h>
int hamming_leq(const char* a, const char* b, int k){
  if (!a || !b) return 0;
  size_t na = strlen(a), nb = strlen(b);
  if (na != nb) return 0;
  int d=0;
  for (size_t i=0;i<na;++i){ if (a[i]!=b[i]){ if (++d>k) return 0; } }
  return 1;
}
