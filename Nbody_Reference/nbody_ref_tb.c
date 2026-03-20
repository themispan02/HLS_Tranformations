//
// Created by themi on 17/3/2026.
//

/*============================================================
 * nbody_naive_tb.c
 *
 * Testbench για τον Naive N-Body Kernel
 *
 * Λειτουργία:
 *   1. Αρχικοποίηση σωματιδίων (τυχαίες θέσεις/ταχύτητες)
 *   2. Εκτέλεση software reference (ίδιος αλγόριθμος)
 *   3. Εκτέλεση HLS kernel
 *   4. Σύγκριση αποτελεσμάτων
 *   5. Έλεγχος διατήρησης ορμής (momentum conservation)
 *
 * Compile (standalone):
 *   gcc -O2 -o nbody_test nbody_naive.c nbody_naive_tb.c -lm
 *
 * Compile (μικρό test):
 *   gcc -O2 -DN_BODIES=128 -DTIMESTEPS=3 -o nbody_test \
 *       nbody_naive.c nbody_naive_tb.c -lm
 *============================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define N_BODIES       16128
#define N_DIMS         3
#define TIMESTEPS      4
#define G_CONST        1.0f
#define SOFTENING_SQ   1000.0f
#define DT             1.0f

/* Πρωτότυπο — ορίζεται στο nbody_naive.c */
void nbody_naive(
    float pos_x[N_BODIES],
    float pos_y[N_BODIES],
    float pos_z[N_BODIES],
    float vel_x[N_BODIES],
    float vel_y[N_BODIES],
    float vel_z[N_BODIES],
    float mass[N_BODIES]
);

/* Ανοχή σύγκρισης — single-precision float */
#define TOLERANCE 1e-3f

/* Πόσα σωματίδια να εκτυπώσουμε για δείγμα */
#define PRINT_SAMPLE 5

/*--- Βοηθητικές Συναρτήσεις ---*/

/* Τυχαίος αριθμός στο [-range, range] */
static float random_range(float range)
{
    return ((float)rand() / (float)RAND_MAX) * 2.0f * range - range;
}

/* Εκτύπωση δείγματος σωματιδίων */
static void print_sample(
    const char *title,
    float *px, float *py, float *pz,
    float *vx, float *vy, float *vz,
    int n, int num_show)
{
    int i;
    printf("\n%s (first %d of %d bodies):\n", title, num_show, n);
    printf("%-6s  %-12s %-12s %-12s  %-12s %-12s %-12s\n",
           "ID", "pos_x", "pos_y", "pos_z", "vel_x", "vel_y", "vel_z");
    printf("------  ------------ ------------ ------------  "
           "------------ ------------ ------------\n");
    for (i = 0; i < num_show && i < n; i++) {
        printf("%-6d  %12.4f %12.4f %12.4f  %12.4f %12.4f %12.4f\n",
               i, px[i], py[i], pz[i], vx[i], vy[i], vz[i]);
    }
    printf("...\n");
}

/* Υπολογισμός συνολικής ορμής */
static void compute_momentum(
    float *vx, float *vy, float *vz,
    float *mass, int n,
    float *mom_x, float *mom_y, float *mom_z)
{
    int i;
    *mom_x = 0.0f; *mom_y = 0.0f; *mom_z = 0.0f;
    for (i = 0; i < n; i++) {
        *mom_x += mass[i] * vx[i];
        *mom_y += mass[i] * vy[i];
        *mom_z += mass[i] * vz[i];
    }
}

/*--- Software Reference (ίδιος αλγόριθμος, timestep-by-timestep) ---*/
static void nbody_reference(
    float *px, float *py, float *pz,
    float *vx, float *vy, float *vz,
    float *mass)
{
    int t, i, j;

    /* Temporary buffers for ping-pong */
    float *npx = (float *)malloc(N_BODIES * sizeof(float));
    float *npy = (float *)malloc(N_BODIES * sizeof(float));
    float *npz = (float *)malloc(N_BODIES * sizeof(float));
    float *nvx = (float *)malloc(N_BODIES * sizeof(float));
    float *nvy = (float *)malloc(N_BODIES * sizeof(float));
    float *nvz = (float *)malloc(N_BODIES * sizeof(float));

    float *spx, *spy, *spz, *svx, *svy, *svz;
    float *dpx, *dpy, *dpz, *dvx, *dvy, *dvz;

    for (t = 0; t < TIMESTEPS; t++) {

        if (t % 2 == 0) {
            spx = px;  spy = py;  spz = pz;
            svx = vx;  svy = vy;  svz = vz;
            dpx = npx; dpy = npy; dpz = npz;
            dvx = nvx; dvy = nvy; dvz = nvz;
        } else {
            spx = npx; spy = npy; spz = npz;
            svx = nvx; svy = nvy; svz = nvz;
            dpx = px;  dpy = py;  dpz = pz;
            dvx = vx;  dvy = vy;  dvz = vz;
        }

        /* Position update */
        for (i = 0; i < N_BODIES; i++) {
            dpx[i] = spx[i] + svx[i] * DT;
            dpy[i] = spy[i] + svy[i] * DT;
            dpz[i] = spz[i] + svz[i] * DT;
        }

        /* Force computation + velocity update */
        for (i = 0; i < N_BODIES; i++) {
            float ax = 0.0f, ay = 0.0f, az = 0.0f;
            float pix = dpx[i], piy = dpy[i], piz = dpz[i];

            for (j = 0; j < N_BODIES; j++) {
                float dx = dpx[j] - pix;
                float dy = dpy[j] - piy;
                float dz = dpz[j] - piz;
                float dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING_SQ;
                float inv_d  = 1.0f / sqrtf(dist_sq);
                float inv_d3 = inv_d * inv_d * inv_d;
                float s = G_CONST * mass[j] * inv_d3;
                ax += s * dx;
                ay += s * dy;
                az += s * dz;
            }

            dvx[i] = svx[i] + ax * DT;
            dvy[i] = svy[i] + ay * DT;
            dvz[i] = svz[i] + az * DT;
        }
    }

    /* Αν TIMESTEPS μονός, αντιγραφή πίσω */
    if (TIMESTEPS % 2 != 0) {
        memcpy(px, npx, N_BODIES * sizeof(float));
        memcpy(py, npy, N_BODIES * sizeof(float));
        memcpy(pz, npz, N_BODIES * sizeof(float));
        memcpy(vx, nvx, N_BODIES * sizeof(float));
        memcpy(vy, nvy, N_BODIES * sizeof(float));
        memcpy(vz, nvz, N_BODIES * sizeof(float));
    }

    free(npx); free(npy); free(npz);
    free(nvx); free(nvy); free(nvz);
}

/*--- Main ---*/
int main()
{
    int i, errors;

    printf("======================================================\n");
    printf("  N-Body Naive Testbench\n");
    printf("  Bodies: %d | Dims: %d | Timesteps: %d\n",
           N_BODIES, N_DIMS, TIMESTEPS);
    printf("  G=%.1f | sf2=%.1f | dt=%.1f\n",
           (float)G_CONST, (float)SOFTENING_SQ, (float)DT);
    printf("======================================================\n");

    /*--- Εκχώρηση μνήμης ---*/
    float *px   = (float *)malloc(N_BODIES * sizeof(float));
    float *py   = (float *)malloc(N_BODIES * sizeof(float));
    float *pz   = (float *)malloc(N_BODIES * sizeof(float));
    float *vx   = (float *)malloc(N_BODIES * sizeof(float));
    float *vy   = (float *)malloc(N_BODIES * sizeof(float));
    float *vz   = (float *)malloc(N_BODIES * sizeof(float));
    float *mass = (float *)malloc(N_BODIES * sizeof(float));

    /* Αντίγραφα για reference */
    float *ref_px = (float *)malloc(N_BODIES * sizeof(float));
    float *ref_py = (float *)malloc(N_BODIES * sizeof(float));
    float *ref_pz = (float *)malloc(N_BODIES * sizeof(float));
    float *ref_vx = (float *)malloc(N_BODIES * sizeof(float));
    float *ref_vy = (float *)malloc(N_BODIES * sizeof(float));
    float *ref_vz = (float *)malloc(N_BODIES * sizeof(float));

    /*--- Αρχικοποίηση: τυχαίες θέσεις/ταχύτητες, μάζα=1 ---*/
    srand(42);
    for (i = 0; i < N_BODIES; i++) {
        px[i] = random_range(100.0f);
        py[i] = random_range(100.0f);
        pz[i] = random_range(100.0f);
        vx[i] = random_range(1.0f);
        vy[i] = random_range(1.0f);
        vz[i] = random_range(1.0f);
        mass[i] = 1.0f;
    }

    /* Αντιγραφή αρχικών τιμών για reference */
    memcpy(ref_px, px, N_BODIES * sizeof(float));
    memcpy(ref_py, py, N_BODIES * sizeof(float));
    memcpy(ref_pz, pz, N_BODIES * sizeof(float));
    memcpy(ref_vx, vx, N_BODIES * sizeof(float));
    memcpy(ref_vy, vy, N_BODIES * sizeof(float));
    memcpy(ref_vz, vz, N_BODIES * sizeof(float));

    /* Αρχική ορμή */
    float init_mx, init_my, init_mz;
    compute_momentum(vx, vy, vz, mass, N_BODIES,
                     &init_mx, &init_my, &init_mz);

    print_sample("Initial state",
                 px, py, pz, vx, vy, vz, N_BODIES, PRINT_SAMPLE);

    /*--- Εκτέλεση Reference ---*/
    printf("\n[1] Running software reference (%d timesteps)...\n",
           TIMESTEPS);
    nbody_reference(ref_px, ref_py, ref_pz,
                    ref_vx, ref_vy, ref_vz, mass);

    print_sample("Reference - Final state",
                 ref_px, ref_py, ref_pz,
                 ref_vx, ref_vy, ref_vz, N_BODIES, PRINT_SAMPLE);

    /*--- Εκτέλεση HLS Kernel ---*/
    printf("\n[2] Running HLS kernel (%d timesteps)...\n",
           TIMESTEPS);
    nbody_naive(px, py, pz, vx, vy, vz, mass);

    print_sample("HLS Kernel - Final state",
                 px, py, pz, vx, vy, vz, N_BODIES, PRINT_SAMPLE);

    /*--- Σύγκριση Αποτελεσμάτων ---*/
    printf("\n[3] Comparing results (tolerance = %.0e)...\n",
           (double)TOLERANCE);
    errors = 0;
    for (i = 0; i < N_BODIES; i++) {
        if (fabsf(px[i] - ref_px[i]) > TOLERANCE ||
            fabsf(py[i] - ref_py[i]) > TOLERANCE ||
            fabsf(pz[i] - ref_pz[i]) > TOLERANCE ||
            fabsf(vx[i] - ref_vx[i]) > TOLERANCE ||
            fabsf(vy[i] - ref_vy[i]) > TOLERANCE ||
            fabsf(vz[i] - ref_vz[i]) > TOLERANCE)
        {
            if (errors < 5) {
                printf("  MISMATCH body %d:\n", i);
                printf("    HLS: pos(%f,%f,%f) vel(%f,%f,%f)\n",
                       px[i], py[i], pz[i], vx[i], vy[i], vz[i]);
                printf("    REF: pos(%f,%f,%f) vel(%f,%f,%f)\n",
                       ref_px[i], ref_py[i], ref_pz[i],
                       ref_vx[i], ref_vy[i], ref_vz[i]);
            }
            errors++;
        }
    }

    /*--- Ελεγχος ορμής ---*/
    float final_mx, final_my, final_mz;
    compute_momentum(vx, vy, vz, mass, N_BODIES,
                     &final_mx, &final_my, &final_mz);

    printf("\n[4] Momentum conservation:\n");
    printf("  Initial: (%e, %e, %e)\n",
           (double)init_mx, (double)init_my, (double)init_mz);
    printf("  Final:   (%e, %e, %e)\n",
           (double)final_mx, (double)final_my, (double)final_mz);
    printf("  Delta:   (%e, %e, %e)\n",
           (double)(final_mx - init_mx),
           (double)(final_my - init_my),
           (double)(final_mz - init_mz));

    /*--- Αποτέλεσμα ---*/
    printf("\n======================================================\n");
    if (errors == 0) {
        printf("  PASS - 0 errors across %d bodies x %d timesteps\n",
               N_BODIES, TIMESTEPS);
    } else {
        printf("  FAIL - %d errors across %d bodies\n",
               errors, N_BODIES);
    }
    printf("======================================================\n");

    /*--- Απελευθέρωση μνήμης ---*/
    free(px);  free(py);  free(pz);
    free(vx);  free(vy);  free(vz);  free(mass);
    free(ref_px); free(ref_py); free(ref_pz);
    free(ref_vx); free(ref_vy); free(ref_vz);

    return (errors == 0) ? 0 : 1;
}