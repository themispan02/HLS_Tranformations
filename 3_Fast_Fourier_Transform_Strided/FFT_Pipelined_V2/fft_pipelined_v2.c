/*
* FFT (strided, radix-2, 1,024 points): pipelined stage.
 *
 * Transformations (§ = section in de Fine Licht et al., TPDS 2021):
 *   §2.6  Loop flattening: the inner loop runs over the N/2 butterflies of a
 *         stage instead of all N indices, and the pair (even, odd) follows from
 *         i through affine expressions. No iteration is wasted and the two
 *         loops are flattened into one pipeline.
 *   §4.1  Memory access extraction: each input is read once into a register
 *         and each output is written once.
 *   DEPENDENCE false: the butterflies of a stage touch disjoint pairs, so the
 *         dependency the tool assumes on the in-place arrays is false.
 *
 * II: 16 (baseline) -> 2, now limited by the ports of the real/img memories.
 */

#include "fft_pipelined_v2.h"

void fft(double real[FFT_SIZE], double img[FFT_SIZE],
         double real_twid[FFT_SIZE/2], double img_twid[FFT_SIZE/2])
{
    outer: for (int log = 0; log < LOG2_SIZE; log++) {

        int span = FFT_SIZE >> (log + 1);

        /* §2.6 One iteration per butterfly */
        inner: for (int i = 0; i < FFT_SIZE / 2; i++) {
            #pragma HLS PIPELINE II=1
            #pragma HLS DEPENDENCE variable=real inter false
            #pragma HLS DEPENDENCE variable=img  inter false

            /* Butterfly pair from i: group = i / span, jj = i % span */
            int group = i >> (LOG2_SIZE - 1 - log);
            int jj    = i &  (span - 1);
            int even  = (group << 1) * span + jj;
            int odd   = even + span;

            /* §4.1 Read once */
            double re = real[even];
            double ro = real[odd];
            double ie = img[even];
            double io = img[odd];

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

            /* §4.1 Write once */
            real[even] = re_new;
            real[odd]  = ro_new;
            img[even]  = ie_new;
            img[odd]   = io_new;
        }
    }
}