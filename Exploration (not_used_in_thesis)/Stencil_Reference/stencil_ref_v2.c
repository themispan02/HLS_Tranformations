
//
// Created by Themis on 15/10/2025.
//

#include <stdio.h>
#include <stdlib.h>

#define N 12  // Διάσταση πίνακα (NxN)
#define ITERATIONS 10  // Αριθμός επαναλήψεων Jacobi

// Τύπος δεδομένων για τους πίνακες
typedef float data_t;

// Top-level function για Vitis HLS - Ενοποιημένη υλοποίηση
void jacobi_stencil(data_t input[N][N], data_t output[N][N]) {
    data_t temp_A[N][N];
    data_t temp_B[N][N];
    int i, j;

    // Αρχικοποίηση του temp_A με τα input δεδομένα
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            temp_A[i][j] = input[i][j];
        }
    }

    // Εκτέλεση επαναλήψεων Jacobi
    for (int iter = 0; iter < ITERATIONS; iter++) {

        // Επιλογή source και destination πινάκων
        data_t (*A)[N] = (iter % 2 == 0) ? temp_A : temp_B;
        data_t (*B)[N] = (iter % 2 == 0) ? temp_B : temp_A;

        // Εφαρμογή του 5-point stencil σε ΟΛΑ τα σημεία
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                // Υπολογισμός γειτόνων με boundary conditions
                const data_t west = (j == 0) ? A[i][j] : A[i][j-1];

                const data_t east = (j == N-1) ? A[i][j] : A[i][j+1];

                const data_t north = (i == 0) ? A[i][j] : A[i-1][j];

                const data_t south = (i == N-1) ? A[i][j] : A[i+1][j];

                // Υπολογισμός νέας τιμής (μέσος όρος 4 γειτόνων)
                B[i][j] = (data_t)0.25 * (north + west + east + south);
            }
        }
    }

    // Αντιγραφή του τελικού αποτελέσματος στο output
    data_t *final = (ITERATIONS % 2 == 0) ? &temp_A[0][0] : &temp_B[0][0];
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            output[i][j] = ((data_t (*)[N])final)[i][j];
        }
    }
}

// Βοηθητική συνάρτηση για εμφάνιση πίνακα
// void print_array(data_t arr[N][N], const char* title) {
//     int i, j;
//     printf("\n%s:\n", title);
//     printf("======================================\n");
//     for (i = 0; i < N; i++) {
//         for (j = 0; j < N; j++) {
//             printf("%6.2f ", arr[i][j]);
//         }
//         printf("\n");
//     }
//     printf("======================================\n\n");
// }
//
// // Testbench για επαλήθευση
// int main() {
//     data_t input[N][N];
//     data_t output[N][N];
//
//     printf("Αρχικοποίηση πίνακα %dx%d...\n", N, N);
//
//     // Αρχικοποίηση με δοκιμαστικά δεδομένα
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < N; j++) {
//             if (i == 0 || i == N-1 || j == 0 || j == N-1) {
//                 // Συνοριακές τιμές (boundary conditions)
//                 input[i][j] = 100.0f;
//             } else {
//                 // Εσωτερικές τιμές
//                 input[i][j] = 0.0f;
//             }
//         }
//     }
//     print_array(input, "ΑΡΧΙΚΟΣ ΠΙΝΑΚΑΣ");
//     printf("Εκτέλεση %d επαναλήψεων Jacobi stencil...\n", ITERATIONS);
//
//     // Κλήση της top-level function
//     jacobi_stencil(input, output, ITERATIONS);
//
//     printf("Ολοκλήρωση!\n\n");
//
//     // Εμφάνιση τελικού πίνακα
//     print_array(output, "ΤΕΛΙΚΟΣ ΠΙΝΑΚΑΣ");
//
//     // Εμφάνιση ενός δείγματος αποτελεσμάτων (γωνίες και κέντρο)
//     printf("Δείγμα αποτελεσμάτων:\n");
//     printf("Πάνω αριστερή γωνία: %.2f\n", output[0][0]);
//     printf("Κέντρο: %.2f\n", output[N/2][N/2]);
//     printf("Κάτω δεξιά γωνία: %.2f\n", output[N-1][N-1]);
//
//     return 0;
// }