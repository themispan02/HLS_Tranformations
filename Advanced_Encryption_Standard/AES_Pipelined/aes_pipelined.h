//
// Created by themi on 3/22/2026.
//

#include <stdint.h>

#define NUM_BLOCKS 64  /* Αριθμός blocks για multi-block ECB */

typedef struct {
    uint8_t key[32];
    uint8_t enckey[32];
    uint8_t deckey[32];
} aes256_context;

/* Single-block (naive, αμετάβλητη) */
void aes256_encrypt_ecb(aes256_context *ctx, uint8_t k[32], uint8_t buf[16]);

////////////////////////////////////////////////////////////////////////////////
// Test harness interface code.
struct bench_args_t {
    aes256_context ctx;
    uint8_t k[32];
    uint8_t buf[NUM_BLOCKS * 16];  /* N blocks × 16 bytes */
};