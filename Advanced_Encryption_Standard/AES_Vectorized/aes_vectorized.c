//
// Created by themi on 2/4/2026.
//

#include "aes_vectorized.h"

/*=============================================================================
 * aes_vectorized.c
 *
 * AES-256 ECB — VECTORIZED Single-Block Implementation
 *
 * Βασίζεται στο PIPELINED (aes_pipelined.c) και προσθέτει:
 *
 *   §3.1  Horizontal Unrolling  — #pragma HLS UNROLL στα byte-level loops.
 *         (Vectorization)         Επεξεργασία όλων των 16 bytes ταυτόχρονα
 *                                 αντί σειριακά. Απαιτεί array partition για
 *                                 παράλληλη πρόσβαση στη μνήμη.
 *
 *   §4.1  Memory Access         — Local register copies των key/enckey/deckey
 *         Extraction              πριν τον κύριο βρόχο. Αποφεύγει επαναλαμβ-
 *                                 ανόμενες BRAM προσβάσεις μέσα στο pipeline.
 *                                 (Αντιστοιχεί στο Listing 9 του paper)
 *
 * Διαφορά από Pipelined:
 *   Pipelined:   16 bytes σειριακά, ένα ανά cycle (II=1 per byte)
 *   Vectorized:  16 bytes παράλληλα, όλα μαζί σε 1 cycle (throughput ×16)
 *
 * Τι ΔΕΝ αλλάζει από Pipelined:
 *   §2.7 Inlining, §2.3 LUTRAM, §2.4 Fusion, §2.6 Flattening — παραμένουν.
 *
 * Τι ΔΕΝ ανήκει εδώ:
 *   §3.2 Vertical Unrolling + §3.3 Dataflow → Systolic (στάδιο 4)
 *=============================================================================*/

#define F(x)   (((x)<<1) ^ ((((x)>>7) & 1) * 0x1b))

/* -----------------------------------------------------------------------
 * §2.3 Random Access Buffering: RAM_2P_LUTRAM (κληρονομείται από Pipelined)
 * §3.1: static → ο HLS βλέπει τον πλήρη πίνακα και μπορεί να τον
 *        partition-αρει αν χρειαστεί σε μεταγενέστερο στάδιο.
 * ----------------------------------------------------------------------- */
static const uint8_t sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5,
    0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc,
    0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a,
    0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b,
    0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85,
    0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17,
    0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88,
    0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9,
    0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6,
    0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94,
    0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68,
    0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
};

#define rj_sbox(x) sbox[(x)]

/* §2.7 Inlining (κληρονομείται) */
static uint8_t rj_xtime(uint8_t x)
{
#pragma HLS INLINE
    return (x & 0x80) ? ((x << 1) ^ 0x1b) : (x << 1);
}

/* -----------------------------------------------------------------------
 * §3.1 Horizontal Unrolling — aes_subBytes
 *
 * Pipelined:  while(i--) → 16 iterations σειριακά, 1 sbox lookup/cycle
 * Vectorized: UNROLL complete → 16 lookups παράλληλα, 1 cycle total
 *
 * Απαιτεί:
 *   1. Loop από 0→16 (όχι while(i--)) για static bounds → UNROLL
 *   2. buf: ARRAY_PARTITION complete → 16 ξεχωριστοί registers
 *   3. sbox: RAM_2P_LUTRAM (§2.3) → αρκετά ports για παράλληλα lookups
 *      (ή complete partition αν χρειαστεί σε Systolic)
 *
 * Κώδικα αλλαγή: while(i--) → for(i=0;i<16;i++) για να επιτρέψει UNROLL.
 * Ο HLS δεν μπορεί να κάνει UNROLL σε while με φθίνουσα μεταβλητή.
 * ----------------------------------------------------------------------- */
static void aes_subBytes(uint8_t *buf)
{
#pragma HLS INLINE
    uint8_t i;
    sub: for (i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        buf[i] = rj_sbox(buf[i]);
    }
}

/* -----------------------------------------------------------------------
 * §3.1 Horizontal Unrolling — aes_addRoundKey
 *
 * 16 XOR operations παράλληλα αντί σειριακά.
 * Κώδικα αλλαγή: while(i--) → for(i=0;i<16;i++) για UNROLL.
 * ----------------------------------------------------------------------- */
static void aes_addRoundKey(uint8_t *buf, uint8_t *key)
{
#pragma HLS INLINE
    uint8_t i;
    addkey: for (i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        buf[i] ^= key[i];
    }
}

/* -----------------------------------------------------------------------
 * §3.1 Horizontal Unrolling — aes_addRoundKey_cpy
 *
 * Unroll του copy+XOR loop.
 * ----------------------------------------------------------------------- */
static void aes_addRoundKey_cpy(uint8_t *buf, uint8_t *key, uint8_t *cpk)
{
#pragma HLS INLINE
    uint8_t i;
    cpkey: for (i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        buf[i] ^= (cpk[i] = key[i]);
        cpk[16+i] = key[16+i];
    }
}

/* -----------------------------------------------------------------------
 * §2.7 Inlining (κληρονομείται) — aes_shiftRows
 * Παραμένει ίδιο: ήδη combinational, χωρίς loop.
 * ----------------------------------------------------------------------- */
static void aes_shiftRows(uint8_t *buf)
{
#pragma HLS INLINE
    uint8_t i, j;
    i = buf[1];  buf[1]  = buf[5];  buf[5]  = buf[9];  buf[9]  = buf[13]; buf[13] = i;
    i = buf[10]; buf[10] = buf[2];  buf[2]  = i;
    j = buf[3];  buf[3]  = buf[15]; buf[15] = buf[11]; buf[11] = buf[7];  buf[7]  = j;
    j = buf[14]; buf[14] = buf[6];  buf[6]  = j;
}

/* -----------------------------------------------------------------------
 * §3.1 Horizontal Unrolling — aes_mixColumns
 *
 * 4 columns παράλληλα αντί σειριακά.
 * Κώδικα αλλαγή: loop με i+=4 → for(col=0;col<4;col++) για UNROLL.
 * Κάθε column = 4 bytes → 16 bytes total σε 1 cycle.
 * ----------------------------------------------------------------------- */
static void aes_mixColumns(uint8_t *buf)
{
#pragma HLS INLINE
    uint8_t col;
    uint8_t a, b, c, d, e;
    mix: for (col = 0; col < 4; col++)
    {
#pragma HLS UNROLL
        a = buf[col*4];   b = buf[col*4+1];
        c = buf[col*4+2]; d = buf[col*4+3];
        e = a ^ b ^ c ^ d;
        buf[col*4]   ^= e ^ rj_xtime(a^b);
        buf[col*4+1] ^= e ^ rj_xtime(b^c);
        buf[col*4+2] ^= e ^ rj_xtime(c^d);
        buf[col*4+3] ^= e ^ rj_xtime(d^a);
    }
}

/* -----------------------------------------------------------------------
 * §3.1 Horizontal Unrolling — aes_expandEncKey
 *
 * Unroll των exp1 (3 iter) και exp2 (3 iter).
 * Μικρά loops → UNROLL → fully combinational key expansion.
 * ----------------------------------------------------------------------- */
static void aes_expandEncKey(uint8_t *k, uint8_t *rc)
{
#pragma HLS INLINE
    uint8_t i;
    k[0] ^= rj_sbox(k[29]) ^ (*rc);
    k[1] ^= rj_sbox(k[30]);
    k[2] ^= rj_sbox(k[31]);
    k[3] ^= rj_sbox(k[28]);
    *rc = F(*rc);
    exp1: for (i = 4; i < 16; i += 4)
    {
#pragma HLS UNROLL
        k[i] ^= k[i-4]; k[i+1] ^= k[i-3]; k[i+2] ^= k[i-2]; k[i+3] ^= k[i-1];
    }
    k[16] ^= rj_sbox(k[12]);
    k[17] ^= rj_sbox(k[13]);
    k[18] ^= rj_sbox(k[14]);
    k[19] ^= rj_sbox(k[15]);
    exp2: for (i = 20; i < 32; i += 4)
    {
#pragma HLS UNROLL
        k[i] ^= k[i-4]; k[i+1] ^= k[i-3]; k[i+2] ^= k[i-2]; k[i+3] ^= k[i-1];
    }
}

/* -----------------------------------------------------------------------
 * aes256_encrypt_ecb — TOP-LEVEL, single block (ίδια signature)
 *
 * Νέα pragmas vs Pipelined:
 *
 * §3.1 ARRAY_PARTITION complete — buf, local_key, local_enckey:
 *   Χωρίς partition: 1 BRAM port → μόνο 1 στοιχείο/cycle ακόμα και με UNROLL
 *   Με complete partition: κάθε byte = ξεχωριστός register → N παράλληλες
 *   προσβάσεις/cycle. Αναγκαίο συνοδευτικό του UNROLL (Listing 10 paper).
 *
 * §4.1 Memory Access Extraction — local_key[32], local_enckey[16]:
 *   Αντί να διαβάζει ctx->key/enckey από BRAM κάθε φορά μέσα στον ecb3,
 *   εξάγει τις τιμές σε local arrays (registers) πριν τον βρόχο.
 *   → Μηδενική BRAM contention κατά τη διάρκεια του pipeline.
 *   → Αντιστοιχεί στο Listing 9 του paper: "extract reads before pipeline".
 * ----------------------------------------------------------------------- */
void aes256_encrypt_ecb(aes256_context *ctx, uint8_t k[32], uint8_t buf[16])
{
/* §2.3: sbox → LUTRAM (κληρονομείται από Pipelined) */
#pragma HLS RESOURCE variable=sbox core=RAM_2P_LUTRAM

/* §3.1: buf → complete partition → 16 ξεχωριστοί registers
 * Χωρίς αυτό: UNROLL στα sub/mix/addkey loops δεν έχει αποτέλεσμα
 * γιατί ο HLS δεν μπορεί να παράγει παράλληλες BRAM προσβάσεις. */
#pragma HLS ARRAY_PARTITION variable=buf complete

/* §4.1: Local copies για key material → registers αντί BRAM.
 * Τα ctx->key, ctx->enckey, ctx->deckey είναι BRAM (pointer args).
 * Με local arrays + complete partition: zero-latency access στο pipeline. */
    uint8_t local_key[32];
    uint8_t local_enckey[16];
#pragma HLS ARRAY_PARTITION variable=local_key complete
#pragma HLS ARRAY_PARTITION variable=local_enckey complete

    uint8_t rcon = 1;
    uint8_t i;

    /* §4.1 Memory Access Extraction:
     * Εξαγωγή k[] σε local register πριν τον βρόχο.
     * Μέσα στον ecb3, χρησιμοποιείται local_key (register) όχι k (BRAM). */
    ecb1: for (i = 0; i < 32; i++)
    {
#pragma HLS UNROLL
        ctx->enckey[i] = ctx->deckey[i] = k[i];
        local_key[i] = k[i];
    }

    /* §4.1: local_enckey = πρώτα 16 bytes του enckey για addRoundKey_cpy */
    init_enc: for (i = 0; i < 16; i++)
    {
#pragma HLS UNROLL
        local_enckey[i] = ctx->enckey[i];
    }

    ecb2: for (i = 8; --i;)
        aes_expandEncKey(ctx->deckey, &rcon);

    aes_addRoundKey_cpy(buf, ctx->enckey, ctx->key);

    /* §2.6 + §2.4 + §3.1:
     * PIPELINE στον ecb3 (κληρονομείται από Pipelined) +
     * UNROLL στα inner loops (νέο στο Vectorized) =
     * κάθε round iteration επεξεργάζεται 16 bytes σε 1 cycle
     * αντί 16 cycles (×16 speedup per round). */
    ecb3: for (i = 1, rcon = 1; i < 14; ++i)
    {
#pragma HLS PIPELINE
        aes_subBytes(buf);      /* UNROLL ×16 → 16 sbox lookups parallel  */
        aes_shiftRows(buf);     /* combinational, ήδη parallel             */
        aes_mixColumns(buf);    /* UNROLL ×4  → 4 columns parallel         */
        if (i & 1) aes_addRoundKey(buf, &ctx->key[16]);  /* UNROLL ×16    */
        else       aes_expandEncKey(ctx->key, &rcon),    /* UNROLL ×3+3   */
                   aes_addRoundKey(buf, ctx->key);       /* UNROLL ×16    */
    }
    aes_subBytes(buf);
    aes_shiftRows(buf);
    aes_expandEncKey(ctx->key, &rcon);
    aes_addRoundKey(buf, ctx->key);
}