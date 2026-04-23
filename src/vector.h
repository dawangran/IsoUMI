#pragma once
#include <stddef.h>
typedef struct { char **data; size_t n, cap; } strvec_t;
void strvec_init(strvec_t *v);
int  strvec_push(strvec_t *v, const char *s);
void strvec_free(strvec_t *v);
