/*
* Molecular dynamics: flattened pipelining variant. Not used in the results.
 *
 * Merges the reset, the neighbour loop and the final reduction into one flat
 * loop over all 4,096 (atom, neighbour) pairs, using conditionals on j.
 * Applies §2.1.2 tiled accumulation (TILE_SIZE 8), §2.6 loop flattening and
 * §4.1 memory access extraction (§ = section in de Fine Licht et al., TPDS 2021).
 *
 * The pipelining stage in the results is MD_Pipelined_V2. The loop structure
 * of this variant is reused inside the processing elements of MD_Systolic.
 */
 
#include "md_pipelined.h"

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
 
    /* §2.6 One flat loop; i and j are advanced by hand, without division */
    int i = 0, j = 0;
 
    flat_loop: for (int iter = 0; iter < nAtoms * maxNeighbors; iter++)
    {
        #pragma HLS PIPELINE II=1
        #pragma HLS DEPENDENCE variable=fx_part inter false
        #pragma HLS DEPENDENCE variable=fy_part inter false
        #pragma HLS DEPENDENCE variable=fz_part inter false
 
        /* §4.1 Position of atom i, read at its first neighbour */
        if (j == 0) {
            i_x = position_x[i];
            i_y = position_y[i];
            i_z = position_z[i];
        }

        /* Neighbour read */
        jidx = NL[i * maxNeighbors + j];
        j_x = position_x[jidx];
        j_y = position_y[jidx];
        j_z = position_z[jidx];
 
        /* Lennard-Jones Force Calculation */
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
 
        /* §2.1.2 The first TILE_SIZE neighbours overwrite their slot, so no reset is needed */
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
 
        /* Last neighbour of atom i: reduce the partial sums and store */
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

        if (j == maxNeighbors - 1) {
            j = 0;
            i++;
        } else {
            j++;
        }
    }
}