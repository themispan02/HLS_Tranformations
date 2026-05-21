//
// Created by themi on 5/8/2026.
//

#include "fft_strided_pipelined.h"

void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    #pragma HLS ARRAY_PARTITION variable=real cyclic factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=img cyclic factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=real_twid cyclic factor=2 dim=1
    #pragma HLS ARRAY_PARTITION variable=img_twid  cyclic factor=2 dim=1

    int even, rootindex;

    outer: for(int stage = 0; stage < LOG2_SIZE; stage++){
        int span = FFT_SIZE >> (stage + 1);

        inner: for(int odd = span; odd < FFT_SIZE; odd++){
            #pragma HLS PIPELINE II=1
            #pragma HLS DEPENDENCE variable=real inter false
            #pragma HLS DEPENDENCE variable=img inter false

            if((odd & span) == 0) continue;

            even = odd ^ span;

            // Local registers — αποφεύγουμε WAR hazard:
            // το twiddle χρειάζεται την αρχική τιμή του real[odd],
            // όχι αυτή που μόλις γράψαμε από το butterfly
            // 1 read έκαστο — αποθηκεύονται σε registers
            double r_e = real[even];
            double r_o = real[odd];
            double i_e = img[even];
            double i_o = img[odd];

            // FIX: υπολογισμός butterfly odd result σε register
            // ΠΡΙ το write — έτσι δεν χρειάζεται επιπλέον read
            double r_on = r_e - r_o;
            double i_on = i_e - i_o;

            // Write even (πάντα χωρίς twiddle)
            real[even] = r_e + r_o;
            img[even]  = i_e + i_o;

            // Twiddle
            rootindex = (even << stage) & (FFT_SIZE - 1);
            if(rootindex){
                double tw_r = real_twid[rootindex];
                double tw_i = img_twid[rootindex];
                double r_on = real[odd];
                double i_on = img[odd];
                real[odd]   = tw_r * r_on - tw_i * i_on;
                img[odd]    = tw_r * i_on + tw_i * r_on;
            }
        }
    }
}