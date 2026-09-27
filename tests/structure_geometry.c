/* Include the implementation to check its private geometric predicate directly. */
#include "../src/structure.c"
#include <assert.h>

int main(void){
  hts_pos_t one_junction[] = {100,102};
  hts_pos_t two_junctions[] = {95,97,105,107};
  observation_t one = {.start=90,.end=120,.junctions=one_junction,.nj=1};
  observation_t two = {.start=90,.end=120,.junctions=two_junctions,.nj=2};
  /* The tolerance overlaps both short introns; a shared chain must still be
     one-to-one. Its explicit-conflict relation must not depend on argument order. */
  assert(!contradicts(&one,&one,10));
  assert(!contradicts(&two,&two,10));
  assert(compatible(&one,&one,10));
  assert(compatible(&two,&two,10));
  assert(contradicts(&one,&two,10));
  assert(contradicts(&two,&one,10));
  assert(!compatible(&one,&two,10));
  assert(!compatible(&two,&one,10));

  hts_pos_t suffix_junctions[] = {105,107,115,117};
  hts_pos_t full_junctions[] = {95,97,105,107,115,117};
  observation_t suffix = {.start=102,.end=125,.junctions=suffix_junctions,.nj=2};
  observation_t full = {.start=90,.end=125,.junctions=full_junctions,.nj=3};
  /* An earlier boundary is also within tolerance, but only offset +1 covers
     every observable junction. A greedy first-match rule would reject this. */
  assert(!contradicts(&suffix,&full,10));
  assert(!contradicts(&full,&suffix,10));
  assert(compatible(&suffix,&full,10));
  assert(compatible(&full,&suffix,10));

  observation_t left = {.start=0,.end=50};
  observation_t right = {.start=100,.end=150};
  assert(!contradicts(&left,&right,10));
  assert(!contradicts(&right,&left,10));
  assert(!compatible(&left,&right,10));
  assert(!compatible(&right,&left,10));

  hts_pos_t intron[] = {100,200};
  observation_t spliced = {.start=50,.end=250,.junctions=intron,.nj=1};
  observation_t partial_retention = {.start=150,.end=250};
  assert(contradicts(&spliced,&partial_retention,10));
  assert(contradicts(&partial_retention,&spliced,10));
  puts("PASS: symmetric structural conflict and nonoverlapping-fragment geometry");
  return 0;
}
