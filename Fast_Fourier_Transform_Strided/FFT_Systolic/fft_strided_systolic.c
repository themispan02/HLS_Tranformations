//
// Created by themi on 6/6/2026.
//

#include "fft_strided_systolic.h"

/* -----------------------------------------------------------------------
 * One butterfly stage.
 * INLINE forces a separate hardware instance per call site, so each of
 * the 10 calls below becomes its own pipeline with a constant log value.
 * With constant log, span is constant -> addresses are affine in i.
 * ----------------------------------------------------------------------- */
static void fft_stage(
    double r_in [FFT_SIZE],       double i_in [FFT_SIZE],
    double r_out[FFT_SIZE],       double i_out[FFT_SIZE],
    double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2],
    int log)
{
    #pragma HLS INLINE
    int span = FFT_SIZE >> (log + 1);

    for (int i = 0; i < FFT_SIZE / 2; i++) {
        #pragma HLS PIPELINE II=1
        /* r_in  is read-only  -> no write conflict                      */
        /* r_out is write-only -> consecutive iterations write to         */
        /* disjoint (even, odd) pairs, so no WAW dependency               */
        #pragma HLS DEPENDENCE variable=r_out inter false
        #pragma HLS DEPENDENCE variable=i_out inter false

        int group = i >> (LOG2_SIZE - 1 - log);
        int jj    = i &  (span - 1);
        int even  = (group << 1) * span + jj;
        int odd   = even + span;

        double re = r_in[even];  double ro = r_in[odd];
        double ie = i_in[even];  double io = i_in[odd];

        double re_new = re + ro;  double ro_new = re - ro;
        double ie_new = ie + io;  double io_new = ie - io;

        int rootindex = (even << log) & (FFT_SIZE - 1);
        if (rootindex) {
            double tr = real_twid[rootindex] * ro_new
                      - img_twid[rootindex]  * io_new;
            io_new    = real_twid[rootindex] * io_new
                      + img_twid[rootindex]  * ro_new;
            ro_new    = tr;
        }

        r_out[even] = re_new;  r_out[odd] = ro_new;
        i_out[even] = ie_new;  i_out[odd] = io_new;
    }
}

/* -----------------------------------------------------------------------
 * Top-level FFT: 10 stages as a DATAFLOW pipeline.
 * Each sNr/sNi buffer has exactly one producer (stage N) and one
 * consumer (stage N+1), satisfying the DATAFLOW canonical form.
 * Stage 9 writes the final result directly back to real[]/img[].
 * ----------------------------------------------------------------------- */
void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    #pragma HLS DATAFLOW

    double s0r[FFT_SIZE], s0i[FFT_SIZE];
    double s1r[FFT_SIZE], s1i[FFT_SIZE];
    double s2r[FFT_SIZE], s2i[FFT_SIZE];
    double s3r[FFT_SIZE], s3i[FFT_SIZE];
    double s4r[FFT_SIZE], s4i[FFT_SIZE];
    double s5r[FFT_SIZE], s5i[FFT_SIZE];
    double s6r[FFT_SIZE], s6i[FFT_SIZE];
    double s7r[FFT_SIZE], s7i[FFT_SIZE];
    double s8r[FFT_SIZE], s8i[FFT_SIZE];

    fft_stage(real, img,  s0r, s0i, real_twid, img_twid, 0);
    fft_stage(s0r,  s0i,  s1r, s1i, real_twid, img_twid, 1);
    fft_stage(s1r,  s1i,  s2r, s2i, real_twid, img_twid, 2);
    fft_stage(s2r,  s2i,  s3r, s3i, real_twid, img_twid, 3);
    fft_stage(s3r,  s3i,  s4r, s4i, real_twid, img_twid, 4);
    fft_stage(s4r,  s4i,  s5r, s5i, real_twid, img_twid, 5);
    fft_stage(s5r,  s5i,  s6r, s6i, real_twid, img_twid, 6);
    fft_stage(s6r,  s6i,  s7r, s7i, real_twid, img_twid, 7);
    fft_stage(s7r,  s7i,  s8r, s8i, real_twid, img_twid, 8);
    fft_stage(s8r,  s8i,  real, img, real_twid, img_twid, 9);
}
