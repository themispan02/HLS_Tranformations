/*
 * 5-point Jacobi stencil, second version: the kernel is merged into
 * jacobi_stencil and the buffers are swapped through pointers.
 * Exploration, not used in the thesis results.
 */

#include <stdio.h>
#include <stdlib.h>

#define N 12
#define ITERATIONS 10

typedef float data_t;

void jacobi_stencil(data_t input[N][N], data_t output[N][N]) {
    data_t temp_A[N][N];
    data_t temp_B[N][N];
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            temp_A[i][j] = input[i][j];
        }
    }

    for (int iter = 0; iter < ITERATIONS; iter++) {

        /* Swap input and output buffers on every iteration */
        data_t (*A)[N] = (iter % 2 == 0) ? temp_A : temp_B;
        data_t (*B)[N] = (iter % 2 == 0) ? temp_B : temp_A;

        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                const data_t west = (j == 0) ? A[i][j] : A[i][j-1];

                const data_t east = (j == N-1) ? A[i][j] : A[i][j+1];

                const data_t north = (i == 0) ? A[i][j] : A[i-1][j];

                const data_t south = (i == N-1) ? A[i][j] : A[i+1][j];

                B[i][j] = (data_t)0.25 * (north + west + east + south);
            }
        }
    }

    data_t *final = (ITERATIONS % 2 == 0) ? &temp_A[0][0] : &temp_B[0][0];
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            output[i][j] = ((data_t (*)[N])final)[i][j];
        }
    }
}
