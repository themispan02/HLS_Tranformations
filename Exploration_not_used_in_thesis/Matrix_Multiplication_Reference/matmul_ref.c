/*
 * Matrix multiplication 64x64 (double): reference. Exploration, not used in
 * the thesis results. The k loop accumulates into a single register, so its
 * floating-point additions form a dependency chain.
 *///
// Created by themi on 12/2/2026. .
//

#include <stdint.h>

#define N 64
#define K 64
#define M 64

void matmul_ref(double A[N][K], double B[K][M], double C[N][M]) {
    for (int n = 0; n < N; n++) {
        for (int m = 0; m < M; m++) {
            double acc = C[n][m];
            #pragma HLS  PIPELINE
            for (int k = 0; k < K; k++) {
                acc += A[n][k] * B[k][m];
            }
            C[n][m] = acc;
        }
    }
}