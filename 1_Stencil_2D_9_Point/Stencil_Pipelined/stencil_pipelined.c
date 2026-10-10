//
// Created by themi on 5/20/2026.
//

#include "stencil_pipelined.h"
/*=============================================================================
 * stencil.c  —  PIPELINED 9-Point Stencil (MachSuite stencil2d)
 *
 * Μετασχηματισμοί από de Fine Licht et al. §2 (Pipeline-Enabling):
 *
 *   §2.2  Delay Buffers
 *         Αντί για 9 reads/cycle από το orig[], κρατάμε 2 line buffers
 *         (delays 1 και 2 rows) + 3×3 sliding window register. Σε κάθε
 *         iteration διαβάζεται ΕΝΑ νέο pixel από το external memory, και
 *         8 από τα 9 στοιχεία του window επαναχρησιμοποιούνται από
 *         προηγούμενες iterations (6 από horizontal shift, 2 από line
 *         buffers). Reduction memory bandwidth: 9×.
 *
 *   §2.6  Pipelined Loop Flattening
 *         Οι δύο εμφωλευμένοι βρόχοι (r, c) μετατρέπονται σε ΕΝΑ flat
 *         loop πάνω από το idx ∈ [0, row_size·col_size). Εξαλείφει το
 *         drain overhead μεταξύ διαδοχικών rows (≈ pipeline_depth ×
 *         row_size wasted cycles στο naive nested approach).
 *
 *   §2.3  Random Access Buffering
 *         Ο filter[9] είναι function argument που κανονικά πάει σε BRAM
 *         με 2 ports. Όταν unrollαρουμε τους k1/k2 βρόχους (trip count 3)
 *         χρειαζόμαστε 9 ταυτόχρονα reads → port contention. Λύση:
 *         αντιγραφή σε local_filter[] με ARRAY_PARTITION complete →
 *         όλα τα στοιχεία σε registers, 9 parallel reads χωρίς conflict.
 *
 * Extras (όχι §2, υποστηρικτικά):
 *   - ARRAY_PARTITION complete στο window/local_filter: απαραίτητο
 *     pragma για να υλοποιηθούν τα παραπάνω σε hardware.
 *   - Tree reduction: αυτόματη από τον HLS αφού οι k1/k2 loops είναι
 *     fully unrolled και TYPE=int32_t (integer add latency = 1).
 *============================================================================*/

void stencil (TYPE orig[row_size * col_size],
              TYPE sol[row_size * col_size],
              TYPE filter[f_size])
{
    /*------------------------------------------------------------------------*
     * §2.2  Delay Buffers — 2 line buffers (BRAM, dual-port)
     *
     * line_buf_0 → πίξελ από row (r-2)  ("top"    γραμμή του window)
     * line_buf_1 → πίξελ από row (r-1)  ("middle" γραμμή του window)
     * Η τρέχουσα γραμμή (r) έρχεται κατευθείαν από το orig[idx].
     *------------------------------------------------------------------------*/
    TYPE line_buf_0[col_size];
    TYPE line_buf_1[col_size];

    /*------------------------------------------------------------------------*
     * §2.2  Delay Buffers — 3×3 sliding window (registers)
     *
     * Σε κάθε iteration κάνουμε horizontal shift: window[*][0] ← window[*][1],
     * window[*][1] ← window[*][2], και το window[*][2] παίρνει τη νέα στήλη
     * (top, mid, bot) από τα line buffers + το new pixel.
     *------------------------------------------------------------------------*/
    TYPE window[3][3];
    #pragma HLS ARRAY_PARTITION variable=window complete dim=0

    /*------------------------------------------------------------------------*
     * §2.3  Random Access Buffering — filter coefficients σε registers
     *------------------------------------------------------------------------*/
    TYPE local_filter[f_size];
    #pragma HLS ARRAY_PARTITION variable=local_filter complete dim=1

    filter_load: for (int i = 0; i < f_size; i++) {
        #pragma HLS PIPELINE II=1
        local_filter[i] = filter[i];
    }

    /*------------------------------------------------------------------------*
     * §2.6  Pipelined Loop Flattening — single pipelined raster-scan loop
     *
     * Αντί για nested  for(r)  for(c)  ... έχουμε ένα flat loop με μήκος
     * row_size × col_size = 8192 iterations. Counters r, c ενημερώνονται
     * manually στο τέλος κάθε iteration.
     *------------------------------------------------------------------------*/
    int r = 0, c = 0;

    main_loop: for (int idx = 0; idx < row_size * col_size; idx++) {
        #pragma HLS PIPELINE II=1

        /*--- Horizontal shift of 3×3 window (data reuse: 6/9 elements) ---*/
        window[0][0] = window[0][1]; window[0][1] = window[0][2];
        window[1][0] = window[1][1]; window[1][1] = window[1][2];
        window[2][0] = window[2][1]; window[2][1] = window[2][2];

        /*--- ONE read from external memory per cycle (vs 9 στο naive) ---*/
        TYPE new_pixel = orig[idx];

        /*--- Νέα δεξιά στήλη του window από τα line buffers ---*/
        TYPE top = line_buf_0[c];      /* πίξελ στο (r-2, c) */
        TYPE mid = line_buf_1[c];      /* πίξελ στο (r-1, c) */
        TYPE bot = new_pixel;          /* πίξελ στο (r,   c) */

        window[0][2] = top;
        window[1][2] = mid;
        window[2][2] = bot;

        /*--- Shift line buffers vertically (top ← mid ← bot) ---*/
        line_buf_0[c] = mid;
        line_buf_1[c] = bot;

        /*--- Emit output όταν το window έχει γεμίσει πλήρως ---*/
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

        /*--- Advance flat counters ---*/
        if (c == col_size - 1) { c = 0; r++; }
        else                   { c++;       }
    }
}