/*
 * Matrix multiplication 64x64 (double): vectorized. Exploration, not used in
 * the thesis results.
 *
 * Adds §3.1 horizontal unrolling to the pipelined version: W elements of the
 * row per iteration, with B and acc reshaped by W so that W values are read in
 * one access. The write-back to C moves to a separate loop.
 */

#include <stdint.h>

#define N 64
#define K 64
#define M 64
#define W 16

void matmul_vectorized(double A[N][K], double B[K][M], double C[N][M]) {
    /* §3.1 W consecutive columns of B in one wide word */
    #pragma HLS ARRAY_RESHAPE variable=B cyclic factor=W dim=2
    
    for (int n = 0; n < N; n++) {
        double acc[M];
        #pragma HLS ARRAY_RESHAPE variable=acc cyclic factor=W dim=1
        
        for (int k = 0; k < K; k++) {
            double a_val = A[n][k];
            
            for (int m = 0; m < M; m += W) {
                #pragma HLS PIPELINE II=1

                /* §3.1 W multiply-accumulates per iteration */
                for (int w = 0; w < W; w++) {
                    #pragma HLS UNROLL
                    
                    int m_idx = m + w;
                    
                    double prev;
                    if (k == 0) {
                        prev = 0.0;
                    } else {
                        prev = acc[m_idx];
                    }
                    
                    double result = prev + a_val * B[k][m_idx];
                    acc[m_idx] = result;
                }
            }
        }
        
        /* Write-back kept out of the compute loop */
        for (int m = 0; m < M; m++) {
            #pragma HLS PIPELINE II=1
            C[n][m] = acc[m];
        }
    }
}