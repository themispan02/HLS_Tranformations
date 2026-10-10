#include "fft_strided_og.h"

void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    outer: for (int log = 0; log < LOG2_SIZE; log++) {

        int span = FFT_SIZE >> (log + 1);   /* 512, 256, 128, ..., 1 */

        inner: for (int i = 0; i < FFT_SIZE / 2; i++) {

            /* Recover butterfly pair from linear index i --------------- */
            int group = i >> (LOG2_SIZE - 1 - log);  /* i / span        */
            int jj    = i &  (span - 1);              /* i % span        */
            int even  = (group << 1) * span + jj;
            int odd   = even + span;

            /* Butterfly on real part ----------------------------------- */
            double tr  = real[even] + real[odd];
            real[odd]  = real[even] - real[odd];
            real[even] = tr;

            /* Butterfly on imaginary part ------------------------------ */
            double ti  = img[even] + img[odd];
            img[odd]   = img[even] - img[odd];
            img[even]  = ti;

            /* Twiddle factor (skip for rootindex == 0, i.e., W^0 = 1) - */
            int rootindex = (even << log) & (FFT_SIZE - 1);
            if (rootindex) {
                tr        = real_twid[rootindex] * real[odd]
                          - img_twid[rootindex]  * img[odd];
                img[odd]  = real_twid[rootindex] * img[odd]
                          + img_twid[rootindex]  * real[odd];
                real[odd] = tr;
            }
        }
    }
}