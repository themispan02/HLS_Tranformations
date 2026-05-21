//
// Created by themi on 10/4/2026.
//

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define TYPE double

// Problem Constants
#define nAtoms        25600
#define maxNeighbors  160
#define MAX_ITERATION 1
// LJ coefficients
#define lj1           1.5
#define lj2           2.0

void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms*maxNeighbors]);
////////////////////////////////////////////////////////////////////////////////
// Test harness interface code.

struct bench_args_t {
    TYPE force_x[nAtoms];
    TYPE force_y[nAtoms];
    TYPE force_z[nAtoms];
    TYPE position_x[nAtoms];
    TYPE position_y[nAtoms];
    TYPE position_z[nAtoms];
    int32_t NL[nAtoms*maxNeighbors];
};