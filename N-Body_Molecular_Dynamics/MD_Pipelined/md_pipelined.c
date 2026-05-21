//
// Created by themi on 10/4/2026.
//

/*=============================================================================
 * md_pipelined.c
 *
 * MD/KNN — PIPELINED Implementation (Flattened)
 *
 * Μετασχηματισμοί από de Fine Licht et al.:
 *
 *   §2.1.2 Tiled Accumulation   — Loop-carried dependency στα fx/fy/fz.
 *          Interleaving           T=8 partial accumulators.
 *                                 Iteration j γράφει στο fx_part[j % T].
 *                                 Πρώτο pass (j<T): overwrite (implicit init).
 *                                 Δεύτερο pass (j≥T): accumulate.
 *                                 Απόσταση T=8 ≥ L_dadd=3 → II=1.
 *
 *   §2.6  Loop Flattening        — Ενοποίηση loop_i × loop_j σε ένα
 *                                  flat pipeline (nAtoms × maxNeighbors
 *                                  = 4096 iterations).
 *                                  Pipeline fill/drain μόνο 1 φορά
 *                                  αντί 256 φορές.
 *                                  Init, compute, collapse ενσωματώνονται
 *                                  με conditionals (if j==0, if j==last).
 *
 *   §2.4  Pipelined Loop Fusion  — Init (j==0) και collapse (j==last)
 *                                  συγχωνεύονται μέσα στο pipeline
 *                                  αντί ξεχωριστών loops.
 *                                  Μηδέν drain cycles μεταξύ atoms.
 *
 *   §4.1  Memory Access         — Εξαγωγή position_x/y/z[i] σε local
 *         Extraction              scalars (i_x, i_y, i_z) στο j==0.
 *
 * Τι ΔΕΝ ανήκει εδώ (→ Vectorized στάδιο):
 *   §3.1 Horizontal Unrolling — UNROLL στο collapse
 *   §3.1 Array Partition      — complete partition στα fx_part arrays
 *   §4.3 Memory Striping      — cyclic partition στα position arrays
 *=============================================================================*/
 
#include "md_pipelined.h"
 
/* T: αριθμός partial accumulators.
 * T ≥ latency_dadd (≈3 για double) → II=1.
 * T=8 διαιρεί ακριβώς maxNeighbors=16. */
#define TILE_SIZE 8
 
void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms*maxNeighbors])
{
    TYPE delx, dely, delz, r2inv;
    TYPE r6inv, potential, force, j_x, j_y, j_z;
    TYPE i_x, i_y, i_z;
    int32_t jidx;
 
    /* §2.1.2 Tiled Accumulation Interleaving */
    TYPE fx_part[TILE_SIZE];
    TYPE fy_part[TILE_SIZE];
    TYPE fz_part[TILE_SIZE];
 
    /* §2.6 Loop Flattening:
     * Ενοποίηση loop_i (nAtoms=256) × loop_j (maxNeighbors=16)
     * σε ένα flat pipeline 4096 iterations.
     * Αποφεύγεται pipeline fill/drain ×256 φορές.
     * Counters i,j διατηρούνται manually αντί division/modulo. */
    int i = 0, j = 0;
 
    flat_loop: for (int iter = 0; iter < nAtoms * maxNeighbors; iter++)
    {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=fx_part inter false
#pragma HLS DEPENDENCE variable=fy_part inter false
#pragma HLS DEPENDENCE variable=fz_part inter false
 
        /* §2.4 Pipelined Loop Fusion + §4.1 Memory Access Extraction:
         * Στην αρχή κάθε atom (j==0), φόρτωσε position σε scalars. */
        if (j == 0) {
            i_x = position_x[i];
            i_y = position_y[i];
            i_z = position_z[i];
        }
 
        /* Ανάγνωση γείτονα */
        jidx = NL[i * maxNeighbors + j];
        j_x = position_x[jidx];
        j_y = position_y[jidx];
        j_z = position_z[jidx];
 
        /* Υπολογισμός δύναμης Lennard-Jones */
        delx = i_x - j_x;
        dely = i_y - j_y;
        delz = i_z - j_z;
        r2inv = 1.0 / (delx * delx + dely * dely + delz * delz);
        r6inv = r2inv * r2inv * r2inv;
        potential = r6inv * (lj1 * r6inv - lj2);
        force = r2inv * potential;
 
        TYPE dx = delx * force;
        TYPE dy = dely * force;
        TYPE dz = delz * force;
 
        /* §2.1.2 Tiled Accumulation Interleaving:
         * j < T:  πρώτο pass → overwrite (implicit init, §2.4)
         * j >= T: δεύτερο pass → accumulate
         * Απόσταση μεταξύ εγγραφών στο ίδιο fx_part[t] = T = 8 ≥ L_dadd */
        int t = j % TILE_SIZE;
        if (j < TILE_SIZE) {
            fx_part[t] = dx;
            fy_part[t] = dy;
            fz_part[t] = dz;
        } else {
            fx_part[t] += dx;
            fy_part[t] += dy;
            fz_part[t] += dz;
        }
 
        /* §2.4 Pipelined Loop Fusion:
         * Στο τέλος κάθε atom (j==last), collapse και store.
         * Collapse inline — χωρίς ξεχωριστό loop.
         * Τρέχει κάθε 16 iterations → δεν επηρεάζει II. */
        if (j == maxNeighbors - 1) {
            TYPE fx = fx_part[0] + fx_part[1] + fx_part[2] + fx_part[3]
                    + fx_part[4] + fx_part[5] + fx_part[6] + fx_part[7];
            TYPE fy = fy_part[0] + fy_part[1] + fy_part[2] + fy_part[3]
                    + fy_part[4] + fy_part[5] + fy_part[6] + fy_part[7];
            TYPE fz = fz_part[0] + fz_part[1] + fz_part[2] + fz_part[3]
                    + fz_part[4] + fz_part[5] + fz_part[6] + fz_part[7];
 
            force_x[i] = fx;
            force_y[i] = fy;
            force_z[i] = fz;
        }
 
        /* Manual counter update (§2.6 — αποφυγή division/modulo) */
        if (j == maxNeighbors - 1) {
            j = 0;
            i++;
        } else {
            j++;
        }
    }
}