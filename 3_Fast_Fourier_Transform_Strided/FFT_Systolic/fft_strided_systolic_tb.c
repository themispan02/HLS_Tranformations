//
// Created by themi on 5/8/2026.
//

#include "fft_strided_systolic.h"
#include <math.h>

int main(){
    double data_x[FFT_SIZE];
    double data_y[FFT_SIZE];
    double img[FFT_SIZE];
    double real[FFT_SIZE];
    int i;

    //set up twiddles...
    double typed;
    int n, N;
    N = FFT_SIZE;

    //Pre-calc twiddles
    for(n=0; n<(N>>1); n++){
        typed = (double)(twoPI*n/N);
        real[n] = cos(typed);
        img[n] = (-1.0)*sin(typed);
    }

    //Init data
    for(i=0; i < FFT_SIZE; i++) {
        data_x[i] = (double)(i);
        data_y[i] = (double)(i);
    }

    /*
    printf("Input:\n");
    for(i=0; i < FFT_SIZE; i++) {
        printf("%f : %f\n", data_x[i], data_y[i]);
    }

    printf("fft:\n");
    */

    fft(data_x, data_y, real, img);

    printf("Output:\n");


    for(i=0; i < FFT_SIZE; i++) {
        printf("%.2f : %.2f\n", data_x[i], data_y[i]);
    }

    return 0;
}