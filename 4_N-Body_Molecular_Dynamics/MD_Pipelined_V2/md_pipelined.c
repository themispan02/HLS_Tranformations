/*
 * Molecular dynamics (k-nearest neighbours): pipelined stage.
 *
 * Transformations (§ = section in de Fine Licht et al., TPDS 2021):
 *   §2.1.2 Tiled accumulation: TILE_SIZE partial sums per force component.
 *          Consecutive neighbours add into different slots, so the latency of
 *          the floating-point adder no longer sets the II.
 *   §2.6   Loop flattening: loop_i and loop_j are flattened by the tool into a
 *          single pipeline of 4,096 iterations.
 *   §4.1   Memory access extraction: the position of atom i is read once per
 *          atom, outside the neighbour loop.
 *
 * II: 9 (baseline) -> 1. The cycle count is still higher than the baseline
 * (4,167 vs 2,419), because the tool had already unrolled the baseline's inner
 * loop on its own. The systolic stage recovers the difference.
 *
 * Compared with the flattened variant in ../MD_Pipelined (TILE_SIZE 8, reset
 * and reduction folded into one loop), this version has the same II and
 * throughput (4,167 vs 4,176 cycles) but needs less hardware: 151 vs 178 DSP,
 * 14,725 vs 23,001 FF and 16,775 vs 23,680 LUT. That is why it is the
 * pipelining stage in the results. The processing elements of MD_Systolic
 * still use the flattened variant.
 */

#include "md_pipelined.h"

/* At least the latency of the double adder, and a divisor of maxNeighbors */
#define TILE_SIZE 4

void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms * maxNeighbors])
{
    /* §2.1.2 Partial sums, one register per slot */
    TYPE fx_part[TILE_SIZE];
    TYPE fy_part[TILE_SIZE];
    TYPE fz_part[TILE_SIZE];
    #pragma HLS ARRAY_PARTITION variable=fx_part complete dim=1
    #pragma HLS ARRAY_PARTITION variable=fy_part complete dim=1
    #pragma HLS ARRAY_PARTITION variable=fz_part complete dim=1

    TYPE i_x, i_y, i_z;

    loop_i: for (int32_t i = 0; i < nAtoms; i++) {

        /* §4.1 Read once per atom */
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
            /* After flattening, the reset of fx_part at an atom boundary looks
             * like a conflict with the last update of the previous atom. The
             * reset always overwrites, so the dependency is false (else II 2). */
            #pragma HLS DEPENDENCE variable=fx_part inter false
            #pragma HLS DEPENDENCE variable=fy_part inter false
            #pragma HLS DEPENDENCE variable=fz_part inter false

            /* §2.1.2 Neighbour j adds into slot j % TILE_SIZE */
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

        /* §2.1.2 Reduce the partial sums once per atom */
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