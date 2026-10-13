/*
* Stencil 2D, 9-point: vectorized stage.
 *
 * Builds on the pipelined stage (§2.2 delay buffering, §2.3 random access
 * buffering, §2.6 loop flattening) and adds (§ = section in de Fine Licht et al.,
 * TPDS 2021):
 *   §3.1  Horizontal unrolling: VEC_WIDTH outputs per cycle. The window widens
 *         to 3 x (VEC_WIDTH + 2), and orig, sol and the line buffers are
 *         partitioned cyclically so that VEC_WIDTH elements can be accessed
 *         in the same cycle.
 *
 * II stays 1; the loop runs 128 x 64/4 = 2,048 times instead of 8,192.
 */

#include "stencil_vectorized.h"

#define VEC_WIDTH 4   /* must divide col_size */

void stencil (TYPE orig[row_size * col_size],
              TYPE sol[row_size * col_size],
              TYPE filter[f_size])
{
    /* §3.1 VEC_WIDTH memory banks, one per lane */
    #pragma HLS ARRAY_PARTITION variable=orig cyclic factor=VEC_WIDTH dim=1
    #pragma HLS ARRAY_PARTITION variable=sol  cyclic factor=VEC_WIDTH dim=1

    /* §2.2 Line buffers for rows r-2 and r-1, banked like orig */
    TYPE line_buf_0[col_size];     /* row (r-2) */
    TYPE line_buf_1[col_size];     /* row (r-1) */
    #pragma HLS ARRAY_PARTITION variable=line_buf_0 cyclic factor=VEC_WIDTH dim=1
    #pragma HLS ARRAY_PARTITION variable=line_buf_1 cyclic factor=VEC_WIDTH dim=1

    /* §2.2 + §3.1 Window widened to 3 x (VEC_WIDTH + 2): lane w uses columns w..w+2 */
    TYPE window[3][VEC_WIDTH + 2];
    #pragma HLS ARRAY_PARTITION variable=window complete dim=0

    /* §2.3 Filter coefficients in registers */
    TYPE local_filter[f_size];
    #pragma HLS ARRAY_PARTITION variable=local_filter complete dim=1

    filter_load: for (int i = 0; i < f_size; i++) {
        #pragma HLS PIPELINE II=1
        local_filter[i] = filter[i];
    }

    /* §2.6 Flat loop, now advancing VEC_WIDTH columns per iteration */
    int r = 0, c = 0;

    main_loop: for (int idx = 0; idx < row_size * (col_size / VEC_WIDTH); idx++) {
        #pragma HLS PIPELINE II=1

        /* Keep the last two columns of the previous step as the new left edge */
        for (int k1 = 0; k1 < 3; k1++) {
            #pragma HLS UNROLL
            window[k1][0] = window[k1][VEC_WIDTH];
            window[k1][1] = window[k1][VEC_WIDTH + 1];
        }

        /* VEC_WIDTH new columns: row r from orig, rows r-2 and r-1 from the line buffers */
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

        for (int w = 0; w < VEC_WIDTH; w++) {
            #pragma HLS UNROLL
            window[0][w + 2] = top_pix[w];
            window[1][w + 2] = mid_pix[w];
            window[2][w + 2] = new_pix[w];
        }

        for (int w = 0; w < VEC_WIDTH; w++) {
            #pragma HLS UNROLL
            line_buf_0[c + w] = mid_pix[w];
            line_buf_1[c + w] = new_pix[w];
        }

        /* §3.1 VEC_WIDTH outputs in parallel; lanes left of column 0 stay idle */
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

        if (c == col_size - VEC_WIDTH) { c = 0; r++;           }
        else                            { c += VEC_WIDTH;       }
    }
}