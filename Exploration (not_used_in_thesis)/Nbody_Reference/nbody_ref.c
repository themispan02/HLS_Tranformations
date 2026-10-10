//
// Created by themi on 17/3/2026.
//

/*============================================================
 * nbody_naive.c
 *
 * Naive N-Body Gravitational Simulation — Kernel
 * Vitis HLS σε C, χωρίς κανένα HLS pragma/βελτιστοποίηση
 *
 * Δομή:
 *   for t = 0 .. TIMESTEPS-1          ← εξωτερικός χρονικός βρόχος
 *     for i = 0 .. N-1                ← position update O(N)
 *       pos_new = pos + vel * dt
 *     for i = 0 .. N-1                ← O(N²) πυρήνας
 *       for j = 0 .. N-1
 *         compute gravity(i,j)
 *       vel_new = vel + acc * dt
 *     swap buffers (ping-pong)
 *
 * Symplectic Euler (kick-drift variant):
 *   Πρώτα ενημερώνεται η θέση με την τρέχουσα ταχύτητα,
 *   μετά υπολογίζεται η δύναμη στις νέες θέσεις,
 *   και τέλος ενημερώνεται η ταχύτητα.
 *
 * Παράμετροι ευθυγραμμισμένοι με:
 *   - Xilinx AI Engine N-Body Tutorial
 *   - de Fine Licht et al. (IEEE TPDS 2021)
 *============================================================*/

#include <math.h>

/*--- Παράμετροι Προσομοίωσης (overridable via -D flags) ---*/
#define N_BODIES       160   /* Αριθμός σωματιδίων                    */
#define N_DIMS         3       /* Διαστάσεις (x, y, z)                  */
#define TIMESTEPS      4       /* Αριθμός χρονικών βημάτων               */

/*--- Φυσικές Σταθερές (Xilinx tutorial conventions) ---*/
#define G_CONST        1.0f    /* Βαρυτική σταθερά (κανονικοποιημένη)   */
#define SOFTENING_SQ   1000.0f /* sf² — αποφυγή singularity             */
#define DT             1.0f    /* Χρονικό βήμα (timestep)               */

/*============================================================
 * nbody_naive — Top-level kernel
 *
 * Εκτελεί TIMESTEPS χρονικά βήματα N-Body simulation.
 * Κάθε timestep:
 *   1. Position update:  x_new = x + vx * dt
 *   2. Force computation: O(N²) all-pairs gravity
 *   3. Velocity update:  vx_new = vx + acc_x * dt
 *
 * Χρησιμοποιεί double-buffering (ping-pong) εσωτερικά:
 *   Ζυγό timestep: pos/vel → buf_pos/buf_vel
 *   Μονό timestep: buf_pos/buf_vel → pos/vel
 *
 * Inputs/Outputs (in-place):
 *   pos_x, pos_y, pos_z  — θέσεις
 *   vel_x, vel_y, vel_z  — ταχύτητες
 *   mass                  — μάζες σωματιδίων
 *============================================================*/
void nbody_naive(
    float pos_x[N_BODIES],
    float pos_y[N_BODIES],
    float pos_z[N_BODIES],
    float vel_x[N_BODIES],
    float vel_y[N_BODIES],
    float vel_z[N_BODIES],
    float mass[N_BODIES])
{
    int i, j;

    /* Εσωτερικοί πίνακες ping-pong (double buffering)
     *
     * Ζυγό timestep: διαβάζουμε από pos/vel, γράφουμε σε buf_pos/buf_vel
     * Μονό timestep: διαβάζουμε από buf_pos/buf_vel, γράφουμε σε pos/vel
     *
     * Στο τέλος, αν TIMESTEPS μονός, αντιγράφουμε πίσω.
     */
    float buf_pos_x[N_BODIES];
    float buf_pos_y[N_BODIES];
    float buf_pos_z[N_BODIES];
    float buf_vel_x[N_BODIES];
    float buf_vel_y[N_BODIES];
    float buf_vel_z[N_BODIES];

    /* Δείκτες σε πηγή (src) και προορισμό (dst)
     * για αποφυγή αντιγραφών — απλά εναλλάσσουμε ρόλους */
    float *src_px, *src_py, *src_pz;
    float *src_vx, *src_vy, *src_vz;
    float *dst_px, *dst_py, *dst_pz;
    float *dst_vx, *dst_vy, *dst_vz;

    /*==========================================================
     * ΚΥΡΙΟΣ ΧΡΟΝΙΚΟΣ ΒΡΟΧΟΣ
     *==========================================================*/
    Loop_timestep: for (int t = 0; t < TIMESTEPS; t++) {

        /* Εναλλαγή ping-pong buffers */
        if (t % 2 == 0) {
            src_px = pos_x;   src_py = pos_y;   src_pz = pos_z;
            src_vx = vel_x;   src_vy = vel_y;   src_vz = vel_z;
            dst_px = buf_pos_x; dst_py = buf_pos_y; dst_pz = buf_pos_z;
            dst_vx = buf_vel_x; dst_vy = buf_vel_y; dst_vz = buf_vel_z;
        } else {
            src_px = buf_pos_x; src_py = buf_pos_y; src_pz = buf_pos_z;
            src_vx = buf_vel_x; src_vy = buf_vel_y; src_vz = buf_vel_z;
            dst_px = pos_x;   dst_py = pos_y;   dst_pz = pos_z;
            dst_vx = vel_x;   dst_vy = vel_y;   dst_vz = vel_z;
        }

        /*======================================================
         * ΒΗΜΑ 1: Ενημέρωση θέσης — O(N)
         *
         * pos_new_i = pos_i + vel_i * dt
         *
         * Χρησιμοποιούμε τις ΤΡΕΧΟΥΣΕΣ ταχύτητες (πριν τον
         * υπολογισμό δυνάμεων αυτού του timestep).
         *======================================================*/
        Loop_position: for (i = 0; i < N_BODIES; i++) {
            dst_px[i] = src_px[i] + src_vx[i] * DT;
            dst_py[i] = src_py[i] + src_vy[i] * DT;
            dst_pz[i] = src_pz[i] + src_vz[i] * DT;
        }

        /*======================================================
         * ΒΗΜΑ 2: Υπολογισμός δυνάμεων + ενημέρωση ταχύτητας
         *
         * Για κάθε σωματίδιο i:
         *   - Υπολογίζουμε τη βαρυτική επιτάχυνση από ΟΛΑ τα j
         *   - Ενημερώνουμε τις ταχύτητες: vel_new = vel + acc * dt
         *
         * Αυτός είναι ο O(N²) πυρήνας:
         *   N_BODIES × N_BODIES = 16128² ≈ 260M force calcs/timestep
         *======================================================*/
        Loop_i: for (i = 0; i < N_BODIES; i++) {

            /* Accumulators επιτάχυνσης — αρχικοποίηση στο μηδέν */
            float acc_x = 0.0f;
            float acc_y = 0.0f;
            float acc_z = 0.0f;

            /* Cache θέση σωματιδίου i σε τοπικές μεταβλητές */
            float pi_x = dst_px[i];
            float pi_y = dst_py[i];
            float pi_z = dst_pz[i];

            /*--------------------------------------------------
             * Εσωτερικός βρόχος: βαρυτική αλληλεπίδραση i ← j
             *
             * F_ij = G * m_j * (r_j - r_i) / |r_j - r_i|³
             *
             * Ο softening factor (sf² = 1000) προστίθεται στο
             * dist² ώστε:
             *  - Όταν i == j → dx=dy=dz=0 → dist²=sf² → acc+=0
             *    (λόγω μηδενικού αριθμητή, η συνεισφορά είναι 0)
             *  - Αποφεύγεται division-by-zero για κοντινά σωματίδια
             *--------------------------------------------------*/
            Loop_j: for (j = 0; j < N_BODIES; j++) {

                /* Διάνυσμα απόστασης r_ij */
                float dx = dst_px[j] - pi_x;
                float dy = dst_py[j] - pi_y;
                float dz = dst_pz[j] - pi_z;

                /* Τετράγωνο απόστασης + softening: ||r||² + sf² */
                float dist_sq = dx * dx + dy * dy + dz * dz + SOFTENING_SQ;

                /* 1 / ||r||³  =  (1/√dist_sq)³ */
                float inv_dist  = 1.0f / sqrtf(dist_sq);
                float inv_dist3 = inv_dist * inv_dist * inv_dist;

                /* Συσσώρευση: acc_i += G * m_j * r_ij / ||r||³ */
                float s = G_CONST * mass[j] * inv_dist3;

                acc_x += s * dx;
                acc_y += s * dy;
                acc_z += s * dz;
            }

            /* Ενημέρωση ταχύτητας: v_{t+1} = v_t + acc * dt */
            dst_vx[i] = src_vx[i] + acc_x * DT;
            dst_vy[i] = src_vy[i] + acc_y * DT;
            dst_vz[i] = src_vz[i] + acc_z * DT;
        }

    } /* end Loop_timestep */

    /*==========================================================
     * Αν TIMESTEPS μονός, τα τελικά αποτελέσματα βρίσκονται
     * στους buf_ πίνακες — αντιγραφή πίσω στα pos/vel arrays
     * ώστε ο caller να βρει τα σωστά αποτελέσματα.
     *==========================================================*/
    if (TIMESTEPS % 2 != 0) {
        for (i = 0; i < N_BODIES; i++) {
            pos_x[i] = buf_pos_x[i];
            pos_y[i] = buf_pos_y[i];
            pos_z[i] = buf_pos_z[i];
            vel_x[i] = buf_vel_x[i];
            vel_y[i] = buf_vel_y[i];
            vel_z[i] = buf_vel_z[i];
        }
    }
}