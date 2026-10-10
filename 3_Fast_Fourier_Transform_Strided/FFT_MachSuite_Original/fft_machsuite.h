/*
 * Header of the original MachSuite FFT (fft/strided).
 * The MachSuite test harness (support.h and struct bench_args_t) is left out,
 * because the kernel is driven by fft_machsuite_tb.c instead.
 */
#include <stdio.h>
#include <stdlib.h>

#define FFT_SIZE 1024
#define twoPI 6.28318530717959

void fft(double real[FFT_SIZE], double img[FFT_SIZE], double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2]);
