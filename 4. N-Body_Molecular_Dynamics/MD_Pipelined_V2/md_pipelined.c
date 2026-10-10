//
// Created by themi on 10/4/2026.
//

/*
 * Implementation based on:
 * A. Danalis, G. Marin, C. McCurdy, J. S. Meredith, P. C. Roth, K. Spafford,
 * V. Tipparaju, and J. S. Vetter.
 * The scalable heterogeneous computing (shoc) benchmark suite.
 * In Proceedings of the 3rd Workshop on General-Purpose Computation on
 * Graphics Processing Units, 2010.
 *
 * Optimizations applied (de Fine Licht et al., IEEE TPDS 2021):
 *
 *   §2.1.2 Tiled Accumulation   — The force accumulators fx, fy, fz
 *          Interleaving           suffer a loop-carried dependency with
 *                                 distance 1 in the flattened loop:
 *                                 each iteration reads the value written
 *                                 by the immediately preceding neighbor
 *                                 computation for the same atom, and the
 *                                 FP-add latency is L_fadd ≈ 4 cycles,
 *                                 forcing II = 4.
 *                                 Fix: TILE_SIZE partial accumulators
 *                                 per dimension. Slot t = j % TILE_SIZE
 *                                 is touched every TILE_SIZE iterations,
 *                                 making the dependency distance TILE_SIZE
 *                                 ≥ L_fadd and allowing II = 1.
 *
 *   §2.6  Loop Flattening /    — loop_i (256) and loop_j (16) are
 *         Coalescing             collapsed into a single 4096-iteration
 *                                pipeline, removing 255 redundant
 *                                fill/drain phases (each ~depth cycles).
 *
 *   §4.1  Memory Access        — position_x/y/z[i] are extracted into
 *         Extraction             scalar registers i_x/y/z at j==0 so
 *                                they are read from the BRAM once per
 *                                atom instead of once per neighbor.
 *
 * TILE_SIZE = 8:
 *   • 8 ≥ L_fadd (double) = 4 cycles        → dependency satisfied
 *   • 8 divides maxNeighbors = 16 exactly   → two updates per slot per atom
 *   • 8 × 3 arrays × 8 bytes = 192 bytes of register state
 */

#include "md_pipelined.h"

/* TILE_SIZE: number of partial accumulators per dimension.
 * Constraint: TILE_SIZE >= L_fadd (double) = 4 cycles.
 * Chosen as 4 — the minimum that satisfies the constraint and divides
 * maxNeighbors = 16 evenly (each slot receives exactly 4 updates per atom). */
#define TILE_SIZE 4

void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms * maxNeighbors])
{
    TYPE fx_part[TILE_SIZE];
    TYPE fy_part[TILE_SIZE];
    TYPE fz_part[TILE_SIZE];
    #pragma HLS ARRAY_PARTITION variable=fx_part complete dim=1
    #pragma HLS ARRAY_PARTITION variable=fy_part complete dim=1
    #pragma HLS ARRAY_PARTITION variable=fz_part complete dim=1

    TYPE i_x, i_y, i_z;

    loop_i: for (int32_t i = 0; i < nAtoms; i++) {

        i_x = position_x[i];
        i_y = position_y[i];
        i_z = position_z[i];

        loop_init: for (int32_t t = 0; t < TILE_SIZE; t++) {
            #pragma HLS UNROLL
            fx_part[t] = (TYPE)0.0;
            fy_part[t] = (TYPE)0.0;
            fz_part[t] = (TYPE)0.0;
        }

        loop_j: for (int32_t j = 0; j < maxNeighbors; j++) {
            #pragma HLS PIPELINE II=1
            /* The tool auto-flattens loop_i and loop_j (triggered by the
             * UNROLL on loop_init), creating a WAW conflict at the atom
             * boundary between the init write and the accumulation write.
             * Both writes target fx_part at distance=1, but the init
             * always unconditionally resets to 0, so the conflict is safe.
             * DEPENDENCE inter false overrides the conservative II=2. */
            #pragma HLS DEPENDENCE variable=fx_part inter false
            #pragma HLS DEPENDENCE variable=fy_part inter false
            #pragma HLS DEPENDENCE variable=fz_part inter false

            int32_t t    = j % TILE_SIZE;
            int32_t jidx = NL[i * maxNeighbors + j];

            TYPE j_x = position_x[jidx];
            TYPE j_y = position_y[jidx];
            TYPE j_z = position_z[jidx];

            TYPE delx = i_x - j_x;
            TYPE dely = i_y - j_y;
            TYPE delz = i_z - j_z;

            TYPE r2inv = (TYPE)1.0 / (delx*delx + dely*dely + delz*delz);
            TYPE r6inv = r2inv * r2inv * r2inv;
            TYPE potential = r6inv * (lj1 * r6inv - lj2);
            TYPE force     = r2inv * potential;

            fx_part[t] += delx * force;
            fy_part[t] += dely * force;
            fz_part[t] += delz * force;
        }

        TYPE fx = (TYPE)0.0;
        TYPE fy = (TYPE)0.0;
        TYPE fz = (TYPE)0.0;
        loop_collapse: for (int32_t t = 0; t < TILE_SIZE; t++) {
            #pragma HLS UNROLL
            fx += fx_part[t];
            fy += fy_part[t];
            fz += fz_part[t];
        }

        force_x[i] = fx;
        force_y[i] = fy;
        force_z[i] = fz;
    }
}