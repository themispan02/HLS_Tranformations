#include <stdint.h>

#define N 64
#define K 64
#define M 64
#define W 16

void matmul_vectorized(double A[N][K], double B[K][M], double C[N][M]) {
    #pragma HLS ARRAY_RESHAPE variable=B cyclic factor=W dim=2
    // ΔΕΝ κάνουμε reshape τον C — δεν τον αγγίζουμε στον υπολογιστικό βρόχο πλέον
    
    for (int n = 0; n < N; n++) {
        double acc[M];
        #pragma HLS ARRAY_RESHAPE variable=acc cyclic factor=W dim=1
        
        // === Κύριος υπολογιστικός βρόχος (μόνο acc και B) ===
        for (int k = 0; k < K; k++) {
            double a_val = A[n][k];
            
            for (int m = 0; m < M; m += W) {
                #pragma HLS PIPELINE II=1
                
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
                    // ΑΦΑΙΡΕΘΗΚΕ: if (k == K-1) { C[n][m_idx] = result; }
                }
            }
        }
        
        // === Ξεχωριστός βρόχος εγγραφής στον C ===
        for (int m = 0; m < M; m++) {
            #pragma HLS PIPELINE II=1
            C[n][m] = acc[m];
        }
    }
}