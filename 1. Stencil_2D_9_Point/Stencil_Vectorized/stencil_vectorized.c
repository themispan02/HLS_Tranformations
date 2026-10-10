//
// Created by themi on 5/20/2026.
//

#include "stencil_vectorized.h"

/*=============================================================================
 * stencil.c  —  VECTORIZED 9-Point Stencil (MachSuite stencil2d)
 *
 * Stage 3 transformations from de Fine Licht et al.:
 *
 *   §3.1  Horizontal Unrolling (Vectorization)
 *         W=4 output pixels υπολογίζονται παράλληλα σε κάθε κύκλο. Το
 *         window διευρύνεται από 3×3 σε 3×(W+2) ώστε να περιέχει όλα τα
 *         pixels που χρειάζονται οι W ταυτόχρονοι MAC trees. Συνολικά
 *         9×W = 36 multipliers ενεργοί ταυτόχρονα.
 *
 * Carried over from Stage 2 (Pipelined):
 *
 *   §2.2  Delay Buffers           (line buffers + sliding window)
 *   §2.6  Pipelined Loop Flattening (τώρα με stride VEC_WIDTH)
 *   §2.3  Random Access Buffering (filter σε registers)
 *
 * Supporting HLS directives:
 *
 *   ARRAY_PARTITION cyclic factor=W on orig, sol, line buffers — απαραίτητο
 *   για W parallel reads/writes per cycle. Χωρίς αυτό, port contention στα
 *   external memory interfaces θα ανέβαζε το II.
 *============================================================================*/


#define VEC_WIDTH 4   /* Output pixels per cycle. Τιμές που πρέπει να
                       * διαιρούν το col_size (=64): 2, 4, 8, 16. */

void stencil (TYPE orig[row_size * col_size],
              TYPE sol[row_size * col_size],
              TYPE filter[f_size])
{
    /*------------------------------------------------------------------------*
     * Partition function arguments — W parallel ports για orig/sol
     *------------------------------------------------------------------------*/
    #pragma HLS ARRAY_PARTITION variable=orig cyclic factor=VEC_WIDTH dim=1
    #pragma HLS ARRAY_PARTITION variable=sol  cyclic factor=VEC_WIDTH dim=1

    /*------------------------------------------------------------------------*
     * §2.2  Delay Buffers — Line buffers (cyclic partition για W ports)
     *------------------------------------------------------------------------*/
    TYPE line_buf_0[col_size];     /* row (r-2) */
    TYPE line_buf_1[col_size];     /* row (r-1) */
    #pragma HLS ARRAY_PARTITION variable=line_buf_0 cyclic factor=VEC_WIDTH dim=1
    #pragma HLS ARRAY_PARTITION variable=line_buf_1 cyclic factor=VEC_WIDTH dim=1

    /*------------------------------------------------------------------------*
     * §2.2  Sliding window — Διευρύνθηκε σε 3 × (W+2) για να φιλοξενεί
     *       τα overlap pixels των W output lanes. Output lane w χρησιμοποιεί
     *       window[*][w], window[*][w+1], window[*][w+2].
     *------------------------------------------------------------------------*/
    TYPE window[3][VEC_WIDTH + 2];
    #pragma HLS ARRAY_PARTITION variable=window complete dim=0

    /*------------------------------------------------------------------------*
     * §2.3  Random Access Buffering — Filter coefficients σε registers
     *------------------------------------------------------------------------*/
    TYPE local_filter[f_size];
    #pragma HLS ARRAY_PARTITION variable=local_filter complete dim=1

    filter_load: for (int i = 0; i < f_size; i++) {
        #pragma HLS PIPELINE II=1
        local_filter[i] = filter[i];
    }

    /*------------------------------------------------------------------------*
     * §2.6 + §3.1: Flat pipelined loop, W outputs per iteration.
     *              Total iterations: row_size × (col_size / VEC_WIDTH)
     *                              = 128 × 16 = 2048   (vs 8192 στο Pipelined)
     *------------------------------------------------------------------------*/
    int r = 0, c = 0;

    main_loop: for (int idx = 0; idx < row_size * (col_size / VEC_WIDTH); idx++) {
        #pragma HLS PIPELINE II=1

        /*--- Shift window LEFT by VEC_WIDTH (κρατάμε 2 leftover στήλες) ---*/
        for (int k1 = 0; k1 < 3; k1++) {
            #pragma HLS UNROLL
            window[k1][0] = window[k1][VEC_WIDTH];
            window[k1][1] = window[k1][VEC_WIDTH + 1];
        }

        /*--- Read W new pixels per row (orig + 2 line buffers) ---*/
        TYPE new_pix[VEC_WIDTH], top_pix[VEC_WIDTH], mid_pix[VEC_WIDTH];
        #pragma HLS ARRAY_PARTITION variable=new_pix complete dim=1
        #pragma HLS ARRAY_PARTITION variable=top_pix complete dim=1
        #pragma HLS ARRAY_PARTITION variable=mid_pix complete dim=1

        for (int w = 0; w < VEC_WIDTH; w++) {
            #pragma HLS UNROLL
            new_pix[w] = orig[r * col_size + c + w];   /* row r           */
            top_pix[w] = line_buf_0[c + w];            /* row (r-2)       */
            mid_pix[w] = line_buf_1[c + w];            /* row (r-1)       */
        }

        /*--- Fill rightmost W columns of window ---*/
        for (int w = 0; w < VEC_WIDTH; w++) {
            #pragma HLS UNROLL
            window[0][w + 2] = top_pix[w];
            window[1][w + 2] = mid_pix[w];
            window[2][w + 2] = new_pix[w];
        }

        /*--- Shift line buffers vertically (top ← mid ← new) ---*/
        for (int w = 0; w < VEC_WIDTH; w++) {
            #pragma HLS UNROLL
            line_buf_0[c + w] = mid_pix[w];
            line_buf_1[c + w] = new_pix[w];
        }

        /*------------------------------------------------------------------*
         * §3.1 Horizontal Unrolling — W outputs computed in parallel
         * Per-lane conditional emit handles το first-cols warm-up cleanly.
         *------------------------------------------------------------------*/
        if (r >= 2) {
            for (int w = 0; w < VEC_WIDTH; w++) {
                #pragma HLS UNROLL
                int out_col = (c - 2) + w;
                if (out_col >= 0) {
                    TYPE temp = 0;
                    stencil_label3: for (int k1 = 0; k1 < 3; k1++) {
                        #pragma HLS UNROLL
                        stencil_label4: for (int k2 = 0; k2 < 3; k2++) {
                            #pragma HLS UNROLL
                            temp += local_filter[k1*3 + k2] * window[k1][w + k2];
                        }
                    }
                    sol[(r - 2) * col_size + out_col] = temp;
                }
            }
        }

        /*--- Advance flat counters: c += W, wrap to next row ---*/
        if (c == col_size - VEC_WIDTH) { c = 0; r++;           }
        else                            { c += VEC_WIDTH;       }
    }
}