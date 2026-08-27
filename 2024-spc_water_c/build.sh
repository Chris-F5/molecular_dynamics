#!/bin/bash

set -ex

CFLAGS="-g -Wall"
LDFLAGS="-lm"

gcc $CFLAGS -c api.c model.c utils.c em.c
gcc $LDFLAGS -o em api.o model.o utils.o em.o
