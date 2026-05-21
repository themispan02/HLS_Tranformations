//
// Created by themi on 5/20/2026.
//

/*=============================================================================
 * stencil_test.c — HLS csim testbench για MachSuite stencil2d 9-point
 *
 * Στρατηγική:
 *   1. Παράγει deterministic input (fixed seed → reproducible runs)
 *   2. Τρέχει software reference (= αυθεντικό MachSuite naive, μετονομασμένο)
 *   3. Τρέχει το DUT (το `stencil` που συνθέτει το HLS)
 *   4. Element-by-element comparison, statistics, PASS/FAIL exit code
 *
 * Το ίδιο testbench χρησιμοποιείται ΑΜΕΤΑΒΛΗΤΟ για κάθε στάδιο
 * (Naive baseline, Pipelined, Vectorized, Systolic). Αλλάζει μόνο το
 * stencil.c που linkάρεται.
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stencil.h"

/*-----------------------------------------------------------------------------
 * Software reference — exact copy του αυθεντικού MachSuite naive kernel.
 * Διαφορετικό όνομα για να συνυπάρχει με το DUT στο ίδιο binary.
 *---------------------------------------------------------------------------*/
static void stencil_ref(TYPE orig  [row_size * col_size],
                        TYPE sol   [row_size * col_size],
                        TYPE filter[f_size])
{
    int r, c, k1, k2;
    TYPE temp, mul;
    for (r = 0; r < row_size - 2; r++) {
        for (c = 0; c < col_size - 2; c++) {
            temp = (TYPE)0;
            for (k1 = 0; k1 < 3; k1++) {
                for (k2 = 0; k2 < 3; k2++) {
                    mul = filter[k1*3 + k2] * orig[(r+k1)*col_size + c+k2];
                    temp += mul;
                }
            }
            sol[r*col_size + c] = temp;
        }
    }
}

int main(void)
{
    /* Buffers — static για να μη φουσκώσει το stack */
    static TYPE orig    [row_size * col_size];
    static TYPE sol_dut [row_size * col_size];
    static TYPE sol_ref [row_size * col_size];
    static TYPE filter  [f_size];

    /*---- Deterministic input generation (fixed seed → reproducible) ----*/
    srand(42);
    for (int i = 0; i < row_size * col_size; i++)
        orig[i] = (TYPE)((rand() % (MAX - MIN + 1)) + MIN);
    for (int i = 0; i < f_size; i++)
        filter[i] = (TYPE)((rand() % (MAX - MIN + 1)) + MIN);

    memset(sol_dut, 0, sizeof(sol_dut));
    memset(sol_ref, 0, sizeof(sol_ref));

    /*---- Run reference (golden) ----*/
    stencil_ref(orig, sol_ref, filter);

    /*---- Run DUT (the kernel HLS will synthesize) ----*/
    stencil(orig, sol_dut, filter);

    /*---- Element-by-element validation στις output θέσεις του naive ----*/
    int       errors        = 0;
    long long sum_abs_diff  = 0;
    TYPE      max_abs_diff  = 0;

    for (int r = 0; r < row_size - 2; r++) {
        for (int c = 0; c < col_size - 2; c++) {
            int  idx      = r * col_size + c;
            TYPE diff     = sol_dut[idx] - sol_ref[idx];
            TYPE abs_diff = (diff < 0) ? -diff : diff;

            if (diff != 0) {
                if (errors < 5) {  /* Show first 5 mismatches */
                    printf("MISMATCH at (r=%d, c=%d): "
                           "DUT=%d  REF=%d  diff=%d\n",
                           r, c,
                           (int)sol_dut[idx], (int)sol_ref[idx], (int)diff);
                }
                errors++;
                sum_abs_diff += abs_diff;
            }
            if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
        }
    }

    /*---- Sanity-check sample (πρώτα 3×3 outputs) ----*/
    printf("\nSample outputs (top-left 3x3):\n");
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            int idx = r * col_size + c;
            printf("  sol[%d,%d] = %d\n", r, c, (int)sol_dut[idx]);
        }
    }

    /*---- Final report ----*/
    int total = (row_size - 2) * (col_size - 2);
    printf("\n========================================\n");
    printf(" Stencil 2D 9-point - Validation Report\n");
    printf("========================================\n");
    printf(" Input size:    %d x %d\n", row_size, col_size);
    printf(" Output size:   %d x %d  (%d elements)\n",
           row_size - 2, col_size - 2, total);
    printf(" Mismatches:    %d / %d  (%.2f%%)\n",
           errors, total, 100.0 * errors / total);
    if (errors > 0) {
        printf(" Max |diff|:    %d\n", (int)max_abs_diff);
        printf(" Mean |diff|:   %.2f\n", (double)sum_abs_diff / errors);
    }
    printf(" RESULT:        %s\n", errors == 0 ? "PASS" : "FAIL");
    printf("========================================\n");

    return (errors == 0) ? 0 : 1;
}