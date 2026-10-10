//
// Created by themi on 5/31/2026.
//
/*=============================================================================
 * md_systolic.c
 *
 * MD/KNN — SYSTOLIC Implementation (P=4 Parallel PEs)
 *
 * Architecture (DATAFLOW με memory-based task handoff, ΧΩΡΙΣ streams):
 *
 *   position_x/y/z ─►[replicate_positions]─► pos_x_p[0..3]
 *                                            pos_y_p[0..3]
 *                                            pos_z_p[0..3]
 *
 *   NL ────────────►[distribute_NL]──────► NL_0, NL_1, NL_2, NL_3
 *
 *   pos_x_p[p], NL_p ─►[md_pe_p] ────────► force_x_p, force_y_p, force_z_p
 *                       (×4 parallel)
 *
 *   force_*_p[0..3] ──►[collect_forces]──► force_x, force_y, force_z
 *
 * Μετασχηματισμοί από de Fine Licht et al.:
 *
 *   §3.2 Vertical Unrolling      — P=4 ανεξάρτητα PEs, καθένα για 64 atoms.
 *                                  Atom-level parallelism (φυσικό για MD γιατί
 *                                  οι atoms είναι αρχιτεκτονικά ανεξάρτητοι).
 *
 *   §3.3 Dataflow               — #pragma HLS DATAFLOW στο top-level.
 *                                  7 ταυτόχρονα tasks: replicate + distribute_NL
 *                                  + 4× md_pe + collect_forces.
 *                                  Memory-based handoff (όχι hls::stream) — ο
 *                                  random access του MD δεν επιτρέπει streaming.
 *
 *   §4.2 Memory Buffering       — Replicated position arrays (4 αντίγραφα).
 *                                  Απαραίτητο γιατί κάθε PE κάνει random access
 *                                  σε ολόκληρα τα positions (μέσω jidx).
 *                                  Cyclic banking δεν εγγυάται conflict-free
 *                                  για random patterns.
 *
 *   §4.3 Memory Striping        — Block partition factor=4 στα NL, force_*.
 *                                  Conflict-free επειδή κάθε PE κάνει block
 *                                  access (sequential range στο NL slice του).
 *
 * Στις μονάδες md_pe διατηρούνται οι μετασχηματισμοί του pipelined:
 *   §2.1.2 Tiled Accumulation Interleaving (T=8 partial accumulators)
 *   §2.6   Loop Flattening (flat 1024-iter loop ανά PE)
 *   §2.4   Pipelined Loop Fusion (init+compute+collapse ενσωματωμένα)
 *   §4.1   Memory Access Extraction (position_x/y/z[global_i] σε scalars)
 *============================================================================*/

#include "md_systolic.h"

/*------------------------------------------------------------------------*
 * TASK 1: replicate_positions
 *
 *   Δημιουργεί P=4 αντίγραφα των position arrays. Απαραίτητο γιατί:
 *   - Κάθε PE κάνει random access pos_x[jidx] με jidx υπολογιζόμενο runtime
 *   - Cyclic partition factor=4 ΔΕΝ εγγυάται conflict-free random access
 *     (4 PEs μπορεί να ζητούν jidx-τιμές που πέφτουν στο ίδιο bank)
 *   - Replication λύνει τελείως το πρόβλημα: κάθε PE βλέπει πλήρες αντίγραφο
 *
 *   Cost: 4 × 256 doubles × 3 arrays = 24 KB (αμελητέο για FPGA)
 *   Latency: 256 cycles + ~5 pipeline depth
 *------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------*
 * TASK 2: distribute_NL
 *
 *   Σπάει το NL[4096] σε 4 blocks των 1024 entries. Κάθε PE βλέπει
 *   το δικό του slice (atoms p*64..(p+1)*64-1 → NL entries p*1024..(p+1)*1024-1).
 *
 *   Με block partition factor=4 στο top-level NL, αυτό γίνεται με 4
 *   ταυτόχρονα reads/cycle από διαφορετικά banks → conflict-free.
 *
 *   Latency: 1024 cycles + pipeline depth
 *------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------*
 * TASKS 3-6: md_pe (καλείται 4 φορές με διαφορετικά atom_start)
 *
 *   Compute kernel για 64 atoms του PE. Πλήρες αντίγραφο του pipelined
 *   με τα ίδια §2.1.2/§2.4/§2.6/§4.1 transforms, αλλά:
 *   - flat_loop trip count = 64 × 16 = 1024 (αντί 4096)
 *   - τοπικοί atom indices [0, 64)
 *   - atom_start argument δίνει την global offset
 *
 *   Διαβάζει: pos_x/y/z (full 256, random access), NL_local (1024, sequential)
 *   Γράφει:  force_x/y/z_local (64, sequential write)
 *
 *   Latency: 1024 cycles + depth (~80) = ~1100 cycles
 *------------------------------------------------------------------------*/
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

    /* §2.1.2 Tiled Accumulation Interleaving */
    TYPE fx_part[TILE_SIZE];
    TYPE fy_part[TILE_SIZE];
    TYPE fz_part[TILE_SIZE];

    /* Manual local counters (§2.6 — αποφυγή division/modulo) */
    int local_i = 0, j = 0;

    flat_loop: for (int iter = 0; iter < ATOMS_PER_PE * maxNeighbors; iter++) {
        #pragma HLS PIPELINE II=1
        #pragma HLS DEPENDENCE variable=fx_part inter false
        #pragma HLS DEPENDENCE variable=fy_part inter false
        #pragma HLS DEPENDENCE variable=fz_part inter false

        /* §2.4 Pipelined Loop Fusion + §4.1 Memory Access Extraction:
         * Στο j==0 (αρχή atom), φόρτωσε position σε scalars.
         * global_i = local_i + atom_start (για σωστό access στις replicated positions) */
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

        /* §2.1.2: implicit init στο πρώτο pass (j<T), accumulate στο δεύτερο */
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

        /* §2.4: Collapse + store inline στο τέλος κάθε atom */
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

        /* Manual counter update */
        if (j == maxNeighbors - 1) {
            j = 0;
            local_i++;
        } else {
            j++;
        }
    }
}

/*------------------------------------------------------------------------*
 * TASK 7: collect_forces
 *
 *   Συναρμολογεί τα 4 partial force arrays (64 entries έκαστο) στα
 *   ενιαία output arrays (256 entries). Block layout:
 *     force_x[  0.. 63] ← in_x_0
 *     force_x[ 64..127] ← in_x_1
 *     force_x[128..191] ← in_x_2
 *     force_x[192..255] ← in_x_3
 *
 *   Με block partition factor=4 στα top-level force outputs, 4 ταυτόχρονα
 *   writes/cycle σε διαφορετικά banks → conflict-free.
 *
 *   Latency: 64 cycles + pipeline depth
 *------------------------------------------------------------------------*/
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

/*------------------------------------------------------------------------*
 * TOP-LEVEL: md_kernel
 *
 *   Instantiates τα 7 tasks σε DATAFLOW region με memory-based handoff.
 *   Τα 4 md_pe instances δημιουργούν ξεχωριστά hardware units (Vitis
 *   παράγει 4 separate instances όταν η ίδια function καλείται σε
 *   DATAFLOW με διαφορετικά arguments).
 *------------------------------------------------------------------------*/
void md_kernel(TYPE force_x[nAtoms],
               TYPE force_y[nAtoms],
               TYPE force_z[nAtoms],
               TYPE position_x[nAtoms],
               TYPE position_y[nAtoms],
               TYPE position_z[nAtoms],
               int32_t NL[nAtoms*maxNeighbors])
{
    /* §4.3 Memory Striping: block partition στα input NL και output force arrays
     * → 4 ports για ταυτόχρονα reads/writes ανά bank, conflict-free */
    #pragma HLS ARRAY_PARTITION variable=NL dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_x dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_y dim=1 type=block factor=4
    #pragma HLS ARRAY_PARTITION variable=force_z dim=1 type=block factor=4

    /* §4.2 Memory Buffering: P=4 αντίγραφα των position arrays
     * complete partition στην dim=1 (P διαφορετικά BRAMs) */
    TYPE pos_x_p[P][nAtoms];
    TYPE pos_y_p[P][nAtoms];
    TYPE pos_z_p[P][nAtoms];
    #pragma HLS ARRAY_PARTITION variable=pos_x_p dim=1 type=complete
    #pragma HLS ARRAY_PARTITION variable=pos_y_p dim=1 type=complete
    #pragma HLS ARRAY_PARTITION variable=pos_z_p dim=1 type=complete

    /* NL slices (4 separate arrays για clean DATAFLOW dataflow) */
    int32_t NL_0[NL_PER_PE], NL_1[NL_PER_PE];
    int32_t NL_2[NL_PER_PE], NL_3[NL_PER_PE];

    /* Force outputs per PE */
    TYPE force_x_0[ATOMS_PER_PE], force_y_0[ATOMS_PER_PE], force_z_0[ATOMS_PER_PE];
    TYPE force_x_1[ATOMS_PER_PE], force_y_1[ATOMS_PER_PE], force_z_1[ATOMS_PER_PE];
    TYPE force_x_2[ATOMS_PER_PE], force_y_2[ATOMS_PER_PE], force_z_2[ATOMS_PER_PE];
    TYPE force_x_3[ATOMS_PER_PE], force_y_3[ATOMS_PER_PE], force_z_3[ATOMS_PER_PE];

    #pragma HLS DATAFLOW

    /* Stage 1: Replicate positions (~260 cycles) */
    replicate_positions(position_x, position_y, position_z,
                        pos_x_p, pos_y_p, pos_z_p);

    /* Stage 2: Distribute NL (~1030 cycles) */
    distribute_NL(NL, NL_0, NL_1, NL_2, NL_3);

    /* Stages 3-6: 4 parallel compute PEs (~1100 cycles each, ταυτόχρονα) */
    md_pe(force_x_0, force_y_0, force_z_0,
          pos_x_p[0], pos_y_p[0], pos_z_p[0], NL_0,   0);
    md_pe(force_x_1, force_y_1, force_z_1,
          pos_x_p[1], pos_y_p[1], pos_z_p[1], NL_1,  64);
    md_pe(force_x_2, force_y_2, force_z_2,
          pos_x_p[2], pos_y_p[2], pos_z_p[2], NL_2, 128);
    md_pe(force_x_3, force_y_3, force_z_3,
          pos_x_p[3], pos_y_p[3], pos_z_p[3], NL_3, 192);

    /* Stage 7: Collect forces (~70 cycles) */
    collect_forces(force_x_0, force_y_0, force_z_0,
                   force_x_1, force_y_1, force_z_1,
                   force_x_2, force_y_2, force_z_2,
                   force_x_3, force_y_3, force_z_3,
                   force_x, force_y, force_z);
}