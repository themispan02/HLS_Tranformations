/*
 * Molecular dynamics (k-nearest neighbours): systolic stage, P = 4 PEs.
 *
 *   position_x/y/z --> replicate_positions --> 4 copies of the positions
 *   NL ------------> distribute_NL -------> one NL slice per PE
 *   copy p, slice p -> md_pe (x4) ---------> forces of 64 atoms each
 *   4 x 64 forces --> collect_forces ------> force_x/y/z
 *
 * Transformations (§ = section in de Fine Licht et al., TPDS 2021):
 *   §3.2  Vertical unrolling: four independent processing elements, each
 *         computing the forces of 64 atoms.
 *   §3.3  Dataflow: the seven tasks above run as concurrent processes,
 *         handing data over through memories (the random neighbour access
 *         rules out streams).
 *   §4.2  Memory buffering: every PE gets its own copy of the positions, since
 *         its neighbour reads can hit any atom.
 *   §4.3  Memory striping: NL and the force outputs are block partitioned by 4,
 *         one bank per PE.
 * Each PE keeps the loop of the flattened MD_Pipelined variant (§2.1.2 tiled
 * accumulation with TILE_SIZE 8, §2.6 loop flattening, §4.1 memory access
 * extraction).
 *
 * Latency 2,198 cycles. The top-level interval is 1,024 cycles, so with a
 * stream of problems a new one starts before the previous one finishes.
 */

#include "md_systolic.h"

/* §4.2 One full copy of the positions per PE */
static void replicate_positions(
    TYPE in_x[nAtoms], TYPE in_y[nAtoms], TYPE in_z[nAtoms],
    TYPE out_x[P][nAtoms], TYPE out_y[P][nAtoms], TYPE out_z[P][nAtoms])
{
    rep_loop: for (int k = 0; k < nAtoms; k++) {
        #pragma HLS PIPELINE II=1
        TYPE vx = in_x[k];
        TYPE vy = in_y[k];
        TYPE vz = in_z[k];

        for (int p = 0; p < P; p++) {
            #pragma HLS UNROLL
            out_x[p][k] = vx;
            out_y[p][k] = vy;
            out_z[p][k] = vz;
        }
    }
}

/* §4.3 PE p gets the neighbour lists of atoms p*64 .. p*64+63, one bank each */
static void distribute_NL(
    int32_t in_NL[nAtoms*maxNeighbors],
    int32_t out_0[NL_PER_PE], int32_t out_1[NL_PER_PE],
    int32_t out_2[NL_PER_PE], int32_t out_3[NL_PER_PE])
{
    dist_loop: for (int k = 0; k < NL_PER_PE; k++) {
        #pragma HLS PIPELINE II=1
        out_0[k] = in_NL[k];
        out_1[k] = in_NL[NL_PER_PE + k];
        out_2[k] = in_NL[2 * NL_PER_PE + k];
        out_3[k] = in_NL[3 * NL_PER_PE + k];
    }
}

/* §3.2 One processing element: 64 atoms, the flattened loop of MD_Pipelined */
static void md_pe(
    TYPE    force_x_local[ATOMS_PER_PE],
    TYPE    force_y_local[ATOMS_PER_PE],
    TYPE    force_z_local[ATOMS_PER_PE],
    TYPE    pos_x[nAtoms],
    TYPE    pos_y[nAtoms],
    TYPE    pos_z[nAtoms],
    int32_t NL_local[NL_PER_PE],
    int     atom_start)
{
    TYPE delx, dely, delz, r2inv;
    TYPE r6inv, potential, force, j_x, j_y, j_z;
    TYPE i_x, i_y, i_z;
    int32_t jidx;

    TYPE fx_part[TILE_SIZE];
    TYPE fy_part[TILE_SIZE];
    TYPE fz_part[TILE_SIZE];

    /* §2.6 One flat loop of 64 x 16 iterations; counters advanced by hand */
    int local_i = 0, j = 0;

    flat_loop: for (int iter = 0; iter < ATOMS_PER_PE * maxNeighbors; iter++) {
        #pragma HLS PIPELINE II=1
        #pragma HLS DEPENDENCE variable=fx_part inter false
        #pragma HLS DEPENDENCE variable=fy_part inter false
        #pragma HLS DEPENDENCE variable=fz_part inter false

        /* §4.1 Position of atom i, read at its first neighbour */
        if (j == 0) {
            int global_i = local_i + atom_start;
            i_x = pos_x[global_i];
            i_y = pos_y[global_i];
            i_z = pos_z[global_i];
        }

        /* Neighbor read */
        jidx = NL_local[local_i * maxNeighbors + j];
        j_x = pos_x[jidx];
        j_y = pos_y[jidx];
        j_z = pos_z[jidx];

        /* Lennard-Jones force calculation */
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

            force_x_local[local_i] = fx;
            force_y_local[local_i] = fy;
            force_z_local[local_i] = fz;
        }

        if (j == maxNeighbors - 1) {
            j = 0;
            local_i++;
        } else {
            j++;
        }
    }
}

/* §4.3 Each PE's 64 forces go to its own bank of the outputs */
static void collect_forces(
    TYPE in_x_0[ATOMS_PER_PE], TYPE in_y_0[ATOMS_PER_PE], TYPE in_z_0[ATOMS_PER_PE],
    TYPE in_x_1[ATOMS_PER_PE], TYPE in_y_1[ATOMS_PER_PE], TYPE in_z_1[ATOMS_PER_PE],
    TYPE in_x_2[ATOMS_PER_PE], TYPE in_y_2[ATOMS_PER_PE], TYPE in_z_2[ATOMS_PER_PE],
    TYPE in_x_3[ATOMS_PER_PE], TYPE in_y_3[ATOMS_PER_PE], TYPE in_z_3[ATOMS_PER_PE],
    TYPE out_x[nAtoms], TYPE out_y[nAtoms], TYPE out_z[nAtoms])
{
    collect_loop: for (int k = 0; k < ATOMS_PER_PE; k++) {
        #pragma HLS PIPELINE II=1
        out_x[k] = in_x_0[k];
        out_y[k] = in_y_0[k];
        out_z[k] = in_z_0[k];

        out_x[ATOMS_PER_PE + k] = in_x_1[k];
        out_y[ATOMS_PER_PE + k] = in_y_1[k];
        out_z[ATOMS_PER_PE + k] = in_z_1[k];

        out_x[2 * ATOMS_PER_PE + k] = in_x_2[k];
        out_y[2 * ATOMS_PER_PE + k] = in_y_2[k];
        out_z[2 * ATOMS_PER_PE + k] = in_z_2[k];

        out_x[3 * ATOMS_PER_PE + k] = in_x_3[k];
        out_y[3 * ATOMS_PER_PE + k] = in_y_3[k];
        out_z[3 * ATOMS_PER_PE + k] = in_z_3[k];
    }
}


void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms*maxNeighbors])
{
    /* §4.3 Block partition by 4: one bank per PE */
    #pragma HLS ARRAY_PARTITION variable=NL dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_x dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_y dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_z dim=1 type=block factor=4

    /* §4.2 Four copies of the positions, each in its own memory */
    TYPE pos_x_p[P][nAtoms];
    TYPE pos_y_p[P][nAtoms];
    TYPE pos_z_p[P][nAtoms];
    #pragma HLS ARRAY_PARTITION variable=pos_x_p dim=1 type=complete
    #pragma HLS ARRAY_PARTITION variable=pos_y_p dim=1 type=complete
    #pragma HLS ARRAY_PARTITION variable=pos_z_p dim=1 type=complete

    int32_t NL_0[NL_PER_PE], NL_1[NL_PER_PE];
    int32_t NL_2[NL_PER_PE], NL_3[NL_PER_PE];

    TYPE force_x_0[ATOMS_PER_PE], force_y_0[ATOMS_PER_PE], force_z_0[ATOMS_PER_PE];
    TYPE force_x_1[ATOMS_PER_PE], force_y_1[ATOMS_PER_PE], force_z_1[ATOMS_PER_PE];
    TYPE force_x_2[ATOMS_PER_PE], force_y_2[ATOMS_PER_PE], force_z_2[ATOMS_PER_PE];
    TYPE force_x_3[ATOMS_PER_PE], force_y_3[ATOMS_PER_PE], force_z_3[ATOMS_PER_PE];

    /* §3.3 Seven concurrent tasks; the four md_pe calls become four hardware instances */
    #pragma HLS DATAFLOW

    replicate_positions(position_x, position_y, position_z,
                        pos_x_p, pos_y_p, pos_z_p);

    distribute_NL(NL, NL_0, NL_1, NL_2, NL_3);

    md_pe(force_x_0, force_y_0, force_z_0,
          pos_x_p[0], pos_y_p[0], pos_z_p[0], NL_0,   0);
    md_pe(force_x_1, force_y_1, force_z_1,
          pos_x_p[1], pos_y_p[1], pos_z_p[1], NL_1,  64);
    md_pe(force_x_2, force_y_2, force_z_2,
          pos_x_p[2], pos_y_p[2], pos_z_p[2], NL_2, 128);
    md_pe(force_x_3, force_y_3, force_z_3,
          pos_x_p[3], pos_y_p[3], pos_z_p[3], NL_3, 192);

    collect_forces(force_x_0, force_y_0, force_z_0,
                   force_x_1, force_y_1, force_z_1,
                   force_x_2, force_y_2, force_z_2,
                   force_x_3, force_y_3, force_z_3,
                   force_x, force_y, force_z);
}