//
// Created by themi on 5/31/2026.
//

/*=============================================================================
 * md_systolic.h  —  MD/KNN Systolic Implementation
 *
 * Interface ίδιο με το pipelined — δεν αλλάζει η function signature.
 * Όλη η αλλαγή είναι εσωτερική (DATAFLOW με P=4 parallel PEs).
 *============================================================================*/

#include <stdint.h>

#define TYPE          double

/* Problem Constants */
#define nAtoms        256
#define maxNeighbors  16
#define MAX_ITERATION 1

/* LJ coefficients */
#define lj1           1.5
#define lj2           2.0

/* Systolic parameters */
#define P             4                              /* Number of parallel PEs */
#define ATOMS_PER_PE  (nAtoms / P)                   /* 64 atoms ανά PE */
#define NL_PER_PE     (ATOMS_PER_PE * maxNeighbors)  /* 1024 NL entries ανά PE */
#define TILE_SIZE     8                              /* §2.1.2 partial accumulators */

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
    TYPE    force_x[nAtoms];
    TYPE    force_y[nAtoms];
    TYPE    force_z[nAtoms];
    TYPE    position_x[nAtoms];
    TYPE    position_y[nAtoms];
    TYPE    position_z[nAtoms];
    int32_t NL[nAtoms*maxNeighbors];
};
