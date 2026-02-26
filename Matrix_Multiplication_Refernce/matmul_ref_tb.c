//
// Created by themi on 26/2/2026.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

#define N 64
#define K 64
#define M 64

// Δήλωση της top-level function
void matmul_ref(double A[N][K], double B[K][M], double C[N][M]);

static double A[N][K];
static double B[K][M];
static double C[N][M];

int main(void) {
    // Αρχικοποίηση
    srand(42);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < K; j++)
            A[i][j] = (double)rand() / RAND_MAX;

    for (int i = 0; i < K; i++)
        for (int j = 0; j < M; j++)
            B[i][j] = (double)rand() / RAND_MAX;

    // Εκτέλεση
    printf("matmul_ref: %dx%dx%d\n", N, K, M);
    clock_t start = clock();
    matmul_ref(A, B, C);
    clock_t end = clock();

    // Αποτελέσματα
    printf("Time: %.6f sec\n", (double)(end - start) / CLOCKS_PER_SEC);
    printf("C[0][0] = %f\n", C[0][0]);
    printf("C[0][1] = %f\n", C[0][1]);
    printf("C[1][0] = %f\n", C[1][0]);
    printf("C[N-1][M-1] = %f\n", C[N-1][M-1]);

    return 0;
}