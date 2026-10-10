/*
 * Original MachSuite FFT (fft/strided), unmodified apart from the header name.
 * Source: https://github.com/breagen/MachSuite/tree/master/fft/strided
 * Copyright (c) Harvard University, BSD-3 license (see MachSuite LICENSE).
 *
 * Kept for reference only. The baseline used in the thesis is ../FFT_Naive,
 * which rewrites the two loop headers so that Vitis HLS can compute the trip
 * counts; the butterfly computation is identical.
 */
#include "fft_machsuite.h"

void fft(double real[FFT_SIZE], double img[FFT_SIZE], double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2]){
    int even, odd, span, log, rootindex;
    double temp;
    log = 0;

    outer:for(span=FFT_SIZE>>1; span; span>>=1, log++){
        inner:for(odd=span; odd<FFT_SIZE; odd++){
            odd |= span;
            even = odd ^ span;

            temp = real[even] + real[odd];
            real[odd] = real[even] - real[odd];
            real[even] = temp;

            temp = img[even] + img[odd];
            img[odd] = img[even] - img[odd];
            img[even] = temp;

            rootindex = (even<<log) & (FFT_SIZE - 1);
            if(rootindex){
                temp = real_twid[rootindex] * real[odd] -
                    img_twid[rootindex]  * img[odd];
                img[odd] = real_twid[rootindex]*img[odd] +
                    img_twid[rootindex]*real[odd];
                real[odd] = temp;
            }
        }
    }
}
