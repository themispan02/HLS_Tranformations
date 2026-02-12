//
// Created by themi on 12/2/2026.
//

#include <stdint.h>

#define N 8192
#define K 8192
#define M 8192

void matmul_naive(double A[N][K], double B[K][M], double C[N][M]) {
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            double acc = C[n][m];
            for (int k = 0; k < K; k++) {
                acc += A[n][k] * B[k][m];
            }
            C[n][m] = acc;
        }
    }
}