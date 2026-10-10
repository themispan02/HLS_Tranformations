 //
// Created by themi on 26/2/2026.
//

#include <stdint.h>

#define N 64
#define K 64
#define M 64

void matmul_pipelined(double A[N][K], double B[K][M], double C[N][M]) {
    for (int n = 0; n < N; n++) {
        // Διάνυσμα (buffer) για τη συσσώρευση μιας ολόκληρης γραμμής.
        double acc[M];

        // Full Transposition 2.1.1: Εναλλαγή των βρόχων k και m
        for (int k = 0; k < K; k++) {
            // Memory Access Extraction 4.1: Διαβάζουμε το A[n][k] μία φορά
            // και το κρατάμε σε καταχωρητή για ολόκληρο τον εσωτερικό βρόχο.
            double a_val = A[n][k];

            for (int m = 0; m < M; m++) {
                // Ενεργοποίηση pipeline με στόχο ένα αποτέλεσμα ανά κύκλο ρολογιού
                #pragma HLS PIPELINE II=1

                // Pipelined Loop Fusion 2.4 & Flattening 2.6: Η αρχικοποίηση, η συσσώρευση και η εγγραφή
                // συγχωνεύονται σε έναν βρόχο για την αποφυγή νεκρών κύκλων.
                double prev;
                if (k == 0) {
                    prev = 0.0; // Αυτόματη αρχικοποίηση στο 0
                } else {
                    prev = acc[m];
                }

                double result = prev + a_val * B[k][m];
                acc[m] = result;
                // Στην τελευταία επανάληψη του k, γράφουμε κατευθείαν το αποτέλεσμα στον πίνακα C
                if (k == K - 1) {
                    C[n][m] = result;
                }
            }
        }
    }
}