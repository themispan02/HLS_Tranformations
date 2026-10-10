#include "fft_strided.h"
#include <math.h>

void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    int even, rootindex;
    double temp;

    // FIX: Γραμμική έξω loop αντί geometric, ώστε HLS να υπολογίσει σωστά τα trip counts
    outer: for(int stage = 0; stage < LOG2_SIZE; stage++){
        int span = FFT_SIZE >> (stage + 1);  // 512, 256, 128, ..., 1
        int log = stage;                       // 0, 1, 2, ..., 9

        inner: for(int odd = span; odd < FFT_SIZE; odd++){
            // Skip iterations που δεν έχουν bit span set
            if((odd & span) == 0) continue;

            even = odd ^ span;

            temp = real[even] + real[odd];
            real[odd] = real[even] - real[odd];
            real[even] = temp;

            temp = img[even] + img[odd];
            img[odd] = img[even] - img[odd];
            img[even] = temp;

            rootindex = (even << log) & (FFT_SIZE - 1);
            if(rootindex){
                temp = real_twid[rootindex] * real[odd] -
                       img_twid[rootindex] * img[odd];
                img[odd] = real_twid[rootindex] * img[odd] +
                           img_twid[rootindex] * real[odd];
                real[odd] = temp;
            }
        }
    }
}