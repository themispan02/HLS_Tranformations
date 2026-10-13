/*
 * 5-point Jacobi stencil on a 12x12 float grid, plain C with a main().
 * Exploration, not used in the thesis results. Borders reuse the centre value.
 */

#include <stdio.h>
#include <stdlib.h>

#define N 12
#define ITERATIONS 10

typedef float data_t;

void jacobi_kernel(data_t A[N][N], data_t B[N][N]) {
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            const data_t west =
                (j == 0) ? A[i][j] : A[i][j-1];

            const data_t east =
                (j == N-1) ? A[i][j] : A[i][j+1];

            const data_t north =
                (i == 0) ? A[i][j] : A[i-1][j];

            const data_t south =
                (i == N-1) ? A[i][j] : A[i+1][j];

            B[i][j] = (data_t)0.25 * (north + west + east + south);
        }
    }
}

void jacobi_stencil(data_t input[N][N], data_t output[N][N], int num_iterations) {
    data_t temp_A[N][N];
    data_t temp_B[N][N];
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            temp_A[i][j] = input[i][j];
        }
    }

    /* Alternate between the two buffers on every iteration */
    for (int iter = 0; iter < num_iterations; iter++) {
        if (iter % 2 == 0) {
            jacobi_kernel(temp_A, temp_B);
        } else {
            jacobi_kernel(temp_B, temp_A);
        }
    }

    data_t *final = (num_iterations % 2 == 0) ? &temp_A[0][0] : &temp_B[0][0];
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            output[i][j] = ((data_t (*)[N])final)[i][j];
        }
    }
}

void print_array(data_t arr[N][N], const char* title) {
    int i, j;
    printf("\n%s:\n", title);
    printf("======================================\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%6.2f ", arr[i][j]);
        }
        printf("\n");
    }
    printf("======================================\n\n");
}

/* Testbench */
int main() {
    data_t input[N][N];
    data_t output[N][N];

    printf("Αρχικοποίηση πίνακα %dx%d...\n", N, N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == 0 || i == N-1 || j == 0 || j == N-1) {
                input[i][j] = 100.0f;
            } else {
                input[i][j] = 0.0f;
            }
        }
    }
    print_array(input, "ΑΡΧΙΚΟΣ ΠΙΝΑΚΑΣ");
    printf("Εκτέλεση %d επαναλήψεων Jacobi stencil...\n", ITERATIONS);

    jacobi_stencil(input, output, ITERATIONS);

    printf("Ολοκλήρωση!\n\n");

    print_array(output, "ΤΕΛΙΚΟΣ ΠΙΝΑΚΑΣ");

    printf("Δείγμα αποτελεσμάτων:\n");
    printf("Πάνω αριστερή γωνία: %.2f\n", output[0][0]);
    printf("Κέντρο: %.2f\n", output[N/2][N/2]);
    printf("Κάτω δεξιά γωνία: %.2f\n", output[N-1][N-1]);

    return 0;
}