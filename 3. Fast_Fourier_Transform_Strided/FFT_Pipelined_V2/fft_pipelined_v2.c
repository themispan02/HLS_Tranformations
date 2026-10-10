//
// Created by themi on 5/8/2026.
//

#include "fft_pipelined_v2.h"

void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    outer: for (int log = 0; log < LOG2_SIZE; log++) {

        int span = FFT_SIZE >> (log + 1);

        inner: for (int i = 0; i < FFT_SIZE / 2; i++) {
            #pragma HLS PIPELINE II=1
            #pragma HLS DEPENDENCE variable=real inter false
            #pragma HLS DEPENDENCE variable=img  inter false

            /* Recover butterfly pair ------------------------------------ */
            int group = i >> (LOG2_SIZE - 1 - log);
            int jj    = i &  (span - 1);
            int even  = (group << 1) * span + jj;
            int odd   = even + span;

            /* Read both inputs once into registers ---------------------- */
            double re = real[even];
            double ro = real[odd];
            double ie = img[even];
            double io = img[odd];

            /* Butterfly (register arithmetic, no extra BRAM access) ----- */
            double re_new = re + ro;
            double ro_new = re - ro;
            double ie_new = ie + io;
            double io_new = ie - io;

            /* Twiddle factor (operates on register values, not re-reads) */
            int rootindex = (even << log) & (FFT_SIZE - 1);
            if (rootindex) {
                double tr = real_twid[rootindex] * ro_new
                          - img_twid[rootindex]  * io_new;
                io_new    = real_twid[rootindex] * io_new
                          + img_twid[rootindex]  * ro_new;
                ro_new    = tr;
            }

            /* Write results once ---------------------------------------- */
            real[even] = re_new;
            real[odd]  = ro_new;
            img[even]  = ie_new;
            img[odd]   = io_new;
        }
    }
}