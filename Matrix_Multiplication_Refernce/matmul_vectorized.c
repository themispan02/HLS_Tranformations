//
// Created by themi on 9/3/2026.
//

#include <stdint.h>

#define N 64
#define K 64
#define M 64

// Ορίζουμε τον παράγοντα διανυσματικοποίησης (Vectorization Factor)
// Θα υπολογίζουμε 16 στοιχεία ανά κύκλο ρολογιού.
#define W 16

void matmul_vectorized(double A[N][K], double B[K][M], double C[N][M]) {
    // Memory Access Transformations (Kastner): 
    // Αυξάνουμε το εύρος της μνήμης ώστε να μπορούμε να διαβάζουμε/γράφουμε W στοιχεία ταυτόχρονα.
    // Εφαρμόζεται στην 2η διάσταση (dim=2) των πινάκων B και C (δηλαδή στις στήλες).
    #pragma HLS ARRAY_RESHAPE variable=B cyclic factor=W dim=2
    #pragma HLS ARRAY_RESHAPE variable=C cyclic factor=W dim=2
    
    for (int n = 0; n < N; n++) {
        double acc[M];
        // Αναδιαμόρφωση και του τοπικού ενταμιευτή (buffer)
        #pragma HLS ARRAY_RESHAPE variable=acc cyclic factor=W dim=1
        
        for (int k = 0; k < K; k++) {
            double a_val = A[n][k];
            
            // Horizontal Unrolling (de Fine Licht 3.1):
            // Strip-mining: Χωρίζουμε τον βρόχο m σε "κομμάτια" (chunks) μεγέθους W
            for (int m = 0; m < M; m += W) {
                // Διατηρούμε τη διοχέτευση (pipeline) στον εξωτερικό strip-mined βρόχο
                #pragma HLS PIPELINE II=1
                
                // Αυτός ο εσωτερικός βρόχος θα "ξετυλιχθεί" πλήρως στο υλικό
                for (int w = 0; w < W; w++) {
                    #pragma HLS UNROLL
                    
                    int m_idx = m + w;
                    
                    // Pipelined Loop Fusion: Η ίδια λογική της if-else για την αρχικοποίηση
                    double prev;
                    if (k == 0) {
                        prev = 0.0;
                    } else {
                        prev = acc[m_idx];
                    }
                    
                    // Οι W πράξεις θα γίνουν ΠΑΡΑΛΛΗΛΑ στον ίδιο κύκλο ρολογιού
                    double result = prev + a_val * B[k][m_idx];
                    acc[m_idx] = result;
                    
                    if (k == K - 1) {
                        C[n][m_idx] = result;
                    }
                }
            }
        }
    }
}