//
// Created by themi on 5/20/2026.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/* Input sizes */
#define col_size 64
#define row_size 128
#define f_size   9

/* Data bounds */
#define TYPE int32_t
#define MAX  1000
#define MIN  1

/* Iterations */
#define MAX_ITERATION 1

void stencil( TYPE orig[row_size * col_size],
              TYPE sol[row_size * col_size],
              TYPE filter[f_size] );

////////////////////////////////////////////////////////////////////////////////
// Test harness interface code.
struct bench_args_t {
    TYPE orig[row_size*col_size];
    TYPE sol[row_size*col_size];
    TYPE filter[f_size];
};
