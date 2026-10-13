/*
* FFT: first pipelining attempt. Not used in the results.
 *
 * Keeps the loop structure of the baseline and adds a pipelined inner loop,
 * DEPENDENCE false and cyclic partitioning by 2. It was replaced by
 * FFT_Pipelined_V2 for two reasons: half of the inner iterations do nothing
 * (the continue), and the odd output is never written back from r_on/i_on, so
 * the result does not match the baseline.
 */

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

            /* Read each input once into a register */
            double r_e = real[even];
            double r_o = real[odd];
            double i_e = img[even];
            double i_o = img[odd];

            double r_on = r_e - r_o;
            double i_on = i_e - i_o;

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