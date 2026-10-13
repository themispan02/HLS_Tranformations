/*
 * Stencil 2D, 9-point: pipelined stage.
 *
 * Transformations (§ = section in de Fine Licht et al., TPDS 2021):
 *   §2.2  Delay buffering: two line buffers and a 3x3 window, so each cycle
 *         reads one new pixel instead of nine.
 *   §2.3  Random access buffering: the filter is copied into registers.
 *   §2.6  Pipelined loop flattening: the row and column loops become a single
 *         loop, so the pipeline does not drain at the end of every row.
 *
 * II: 5 (baseline) -> 1.
 */

#include "stencil_pipelined.h"

void stencil (TYPE orig[row_size * col_size],
              TYPE sol[row_size * col_size],
              TYPE filter[f_size])
{
    /* §2.2 Line buffers: line_buf_0 holds row r-2, line_buf_1 holds row r-1 */
    TYPE line_buf_0[col_size];
    TYPE line_buf_1[col_size];

    /* §2.2 3x3 sliding window, fully partitioned into registers */
    TYPE window[3][3];
    #pragma HLS ARRAY_PARTITION variable=window complete dim=0

    /* §2.3 Filter coefficients in registers, so all nine can be read in one cycle */
    TYPE local_filter[f_size];
    #pragma HLS ARRAY_PARTITION variable=local_filter complete dim=1

    filter_load: for (int i = 0; i < f_size; i++) {
        #pragma HLS PIPELINE II=1
        local_filter[i] = filter[i];
    }

    /* §2.6 One flat loop over the whole grid; r and c are advanced by hand */
    int r = 0, c = 0;

    main_loop: for (int idx = 0; idx < row_size * col_size; idx++) {
        #pragma HLS PIPELINE II=1

        /*--- Horizontal shift of 3×3 window (data reuse: 6/9 elements) ---*/
        window[0][0] = window[0][1]; window[0][1] = window[0][2];
        window[1][0] = window[1][1]; window[1][1] = window[1][2];
        window[2][0] = window[2][1]; window[2][1] = window[2][2];

        /*--- ONE read from external memory per cycle (vs 9 στο naive) ---*/
        TYPE new_pixel = orig[idx];

        /* New right column: rows r-2 and r-1 from the line buffers, row r from the input */
        TYPE top = line_buf_0[c];      /* πίξελ στο (r-2, c) */
        TYPE mid = line_buf_1[c];      /* πίξελ στο (r-1, c) */
        TYPE bot = new_pixel;          /* πίξελ στο (r,   c) */

        window[0][2] = top;
        window[1][2] = mid;
        window[2][2] = bot;

        /* Move the column one row down in the line buffers (top ← mid ← bot) ---*/
        line_buf_0[c] = mid;
        line_buf_1[c] = bot;

        /* Produce an output once the window holds a full 3x3 neighbourhood */
        if (r >= 2 && c >= 2) {
            TYPE temp = 0;
            stencil_label3: for (int k1 = 0; k1 < 3; k1++) {
                #pragma HLS UNROLL
                stencil_label4: for (int k2 = 0; k2 < 3; k2++) {
                    #pragma HLS UNROLL
                    temp += local_filter[k1*3 + k2] * window[k1][k2];
                }
            }
            sol[(r - 2) * col_size + (c - 2)] = temp;
        }

        if (c == col_size - 1) { c = 0; r++; }
        else                   { c++;       }
    }
}