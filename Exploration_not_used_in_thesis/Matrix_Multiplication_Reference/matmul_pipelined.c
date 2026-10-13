/*
 * Matrix multiplication 64x64 (double): pipelined. Exploration, not used in
 * the thesis results.
 *
 * Transformations (§ = section in de Fine Licht et al., TPDS 2021):
 *   §2.1.1 Full transposition: the k and m loops are swapped, so consecutive
 *          iterations accumulate into different elements of acc.
 *   §4.1   Memory access extraction: A[n][k] is read once per k.
 *   §2.4   Pipelined loop fusion: initialisation and write-back are folded
 *          into the main loop with conditionals on k.
 */

#include <stdint.h>

#define N 64
#define K 64
#define M 64

void matmul_pipelined(double A[N][K], double B[K][M], double C[N][M]) {
    for (int n = 0; n < N; n++) {
        /* One accumulator per element of row n */
        double acc[M];

        /* §2.1.1 k outside, m inside */
        for (int k = 0; k < K; k++) {
            /* §4.1 */
            double a_val = A[n][k];

            for (int m = 0; m < M; m++) {
                // Ενεργοποίηση pipeline με στόχο ένα αποτέλεσμα ανά κύκλο ρολογιού
                #pragma HLS PIPELINE II=1

                /* §2.4 Reset at k == 0, write back at k == K-1 */
                double prev;
                if (k == 0) {
                    prev = 0.0; // Αυτόματη αρχικοποίηση στο 0
                } else {
                    prev = acc[m];
                }

                double result = prev + a_val * B[k][m];
                acc[m] = result;
                if (k == K - 1) {
                    C[n][m] = result;
                }
            }
        }
    }
}