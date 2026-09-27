#pragma once
#include <stddef.h>
#include <htslib/sam.h>

/* Temporary per-bucket assignments. Ordinals include passthrough records. */
typedef struct {
  void* observations;
  size_t n, cap;
  char** keys;
  unsigned char* statuses;
  size_t ordinal_count;
  char** owned_keys;
  size_t group_count;
} structure_map_t;

/* scope is transferred; neighborhood/raw remain borrowed until finish. */
int structure_map_add(structure_map_t* map, char* scope,
                      const char* neighborhood, const char* raw,
                      const bam1_t* read, long ordinal);
int structure_map_finish(structure_map_t* map, int tolerance, int min_support);
const char* structure_map_key(const structure_map_t* map, long ordinal);
const char* structure_map_status(const structure_map_t* map, long ordinal);
void structure_map_destroy(structure_map_t* map);
