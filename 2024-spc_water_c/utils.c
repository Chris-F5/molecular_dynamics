#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "utils.h"

void *
xmalloc(size_t size)
{
  void *p;
  if (!(p = malloc(size))) {
    fprintf(stderr, "malloc failed: %s\n", strerror(errno));
    exit(1);
  }
  return p;
}
