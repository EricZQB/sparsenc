/************************************************************************
 * galois_neon_32.c
 * Functions of Galois field arithmetic.
 ************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "arm_neon.h"
#include "galois.h" 
#define GF_POWER    8
static int constructed = 0;
static uint8_t galois_log_table[1<<GF_POWER];
static uint8_t galois_ilog_table[(1<<GF_POWER)];
static uint8_t galois_mult_table[(1<<GF_POWER)*(1<<GF_POWER)];
static uint8_t galois_divi_table[(1<<GF_POWER)*(1<<GF_POWER)];

/* Two half tables are used for NEON multiply_add_region*/
#if defined(ARM_NEON32) || defined(ARM_NEON64)
static uint8_t galois_half_mult_table_high[(1<<GF_POWER)][(1<<(GF_POWER/2))];
static uint8_t galois_half_mult_table_low[(1<<GF_POWER)][(1<<(GF_POWER/2))];
#endif

static int primitive_poly_8  = 0435;    /* 100 011 101: x^8 + x^4 + x^3 + x^2 + 1 */
static int galois_create_log_table();
static int galois_create_mult_table();

int GFConstructed() {
    return constructed;
}

int constructField()
{
    if (constructed)
        return 0;
    else {
        if (galois_create_mult_table() < 0) {
            perror("constructField");
            exit(1);
        }
#if defined(ARM_NEON32) || defined(ARM_NEON64)
    /*
     * Create half tables for NEON multiply_add_region:
     * low table contains the products of an element with all 4-bit words;
     * high table contains the products of an element with all 8-bit words
     * whose last 4 bits are all zero. So each half table contains 256 rows
     * and 16 columns.
     */
    int a, b, c, d;
    int pp = primitive_poly_8;
    for (a = 1; a < (1<<(GF_POWER/2)) ; a++) {
        b = 1;
        c = a;
        d = (a << (GF_POWER/2));
        do {
            galois_half_mult_table_low[b][a] = c;
            galois_half_mult_table_high[b][a] = d;
            b <<= 1;
            if (b & (1<<GF_POWER)) b ^= pp;
            c <<= 1;
            if (c & (1<<GF_POWER)) c ^= pp;
            d <<= 1;
            if (d & (1<<GF_POWER)) d ^= pp;
        } while (c != a);
    }
#endif
        constructed = 1;
    }
    return 0;
}

static int galois_create_log_table()
{
    int j, b;
    int m = GF_POWER;

    int gf_poly = primitive_poly_8;
    int nw      =  1 << GF_POWER;
    int nwml    = (1 << GF_POWER) - 1;

    for (j=0; j<nw; j++) {
        galois_log_table[j] = nwml;
        galois_ilog_table[j] = 0;
    }

    b = 1;
    for (j=0; j<nwml; j++) {
        if (galois_log_table[b] != nwml) {
            fprintf(stderr, "Galois_create_log_tables Error: j=%d, b=%d, B->J[b]=%d, J->B[j]=%d (0%o)\n", j, b, galois_log_table[b], galois_ilog_table[j], (b << 1) ^ gf_poly);
            exit(1);
        }
        galois_log_table[b] = j;
        galois_ilog_table[j] = b;
        b = b << 1;
        if (b & nw)
            b = (b ^ gf_poly) & nwml;
    }

    return 0;
}

static int galois_create_mult_table()
{
    int j, x, y, logx;
    int nw = (1<<GF_POWER);

    // create tables
    if (galois_create_log_table() < 0) {
        fprintf(stderr, "create log/ilog tables failed\n");
        return -1;
    }

    /* Set mult/div tables for x = 0 */
    j = 0;
    galois_mult_table[j] = 0;   /* y = 0 */
    galois_divi_table[j] = -1;
    j++;
    for (y=1; y<nw; y++) {   /* y > 0 */
        galois_mult_table[j] = 0;
        galois_divi_table[j] = 0;
        j++;
    }

    for (x=1; x<nw; x++) {  /* x > 0 */
        galois_mult_table[j] = 0; /* y = 0 */
        galois_divi_table[j] = -1;
        j++;
        logx = galois_log_table[x];

        for (y=1; y<nw; y++) {  /* y > 0 */
            int tmp;
            tmp = logx + galois_log_table[y];
            if (tmp >= ((1<<GF_POWER) - 1))
                tmp -= ((1<<GF_POWER) - 1);             // avoid cross the boundary of log/ilog tables
            galois_mult_table[j] = galois_ilog_table[tmp];

            tmp = logx - galois_log_table[y];
            while (tmp < 0)
                tmp += ((1<<GF_POWER) - 1);
            galois_divi_table[j] = galois_ilog_table[tmp];

            j++;
        }
    }

    return 0;
}

// add operation over GF(2^m)
inline uint8_t galois_add(uint8_t a, uint8_t b)
{
    return a ^ b;
}

inline uint8_t galois_sub(uint8_t a, uint8_t b)
{
    return a ^ b;
}

inline uint8_t galois_multiply(uint8_t a, uint8_t b)
{
    if (a ==0 || b== 0)
        return 0;

    if (a == 1)
        return b;
    else if (b == 1)
        return a;

    uint8_t result = galois_mult_table[(a<<GF_POWER) | b];
    return result;
}

// return a/b
inline uint8_t galois_divide(uint8_t a, uint8_t b)
{
    if (b == 0) {
        fprintf(stderr, "ERROR! Divide by ZERO!\n");
        return -1;
    }

    if (a == 0)
        return 0;

    if (b == 1)
        return a;

    uint8_t result =  galois_divi_table[(a<<GF_POWER) | b];
    return result;
}

/*
 * When NEON is enabled, use NEON instructions to do multiply_add_region
 */
void galois_multiply_add_region(uint8_t *dst, uint8_t *src, uint8_t multiplier, int bytes)
{
    if (multiplier == 0) {
        // add nothing to bytes starting from *dst, just return
        return;
    }
    int i;
#if defined(ARM_NEON32)
    uint8_t *sptr, *dptr, *top;
    sptr = src;
    dptr = dst;
    top  = src + bytes;

    uint8_t *bh, *bl;
    uint8x16_t loset;
    uint8x8x2_t mth, mtl;

    if (multiplier != 1) {
        /* half tables only needed for multiplier != 1 */
        bh = (uint8_t*) galois_half_mult_table_high;
        bh += (multiplier << 4);
        bl = (uint8_t*) galois_half_mult_table_low;
        bl += (multiplier << 4);
        // read split tables as 128-bit values
        mth.val[0] = vld1_u8(bh);
        mth.val[1] = vld1_u8(bh+8);
        mtl.val[0] = vld1_u8(bl);
        mtl.val[1] = vld1_u8(bl+8);
		    loset = vdupq_n_u8(0x0f);
	}

    uint8x16_t va, vb, r, r2, t1;
    uint8x8_t sh, sl, rh, rl;
    while (sptr < top)
    {

        if (sptr + 16 > top) {
            /* remaining data doesn't fit into __m128i, do not use SSE */
            for (i=0; i<top-sptr; i++) {
                if (multiplier == 1)
                    *(dptr+i) ^= *(sptr+i);
                else
                    *(dptr+i) ^= galois_mult_table[((*(sptr+i))<<GF_POWER) | multiplier];
            }
            break;
        }


        va = vld1q_u8 (sptr);
        if (multiplier == 1) {
            /* just XOR */
            vb = vld1q_u8 (dptr);
            vb = veorq_u8(va, vb);
            vst1q_u8 (dptr, vb);
        } else {
            /* use half tables */
            t1 = vandq_u8 (loset, va);    // obtain lower 4-bit of the 16 src elements
            sh = vget_high_u8(t1);
            sl = vget_low_u8(t1);
            rh = vtbl2_u8 (mtl, sh);    // obtain products of the lower 4-bit 
            rl = vtbl2_u8 (mtl, sl);
            r  = vcombine_u8 (rl, rh);
            
            va = vshrq_n_u8 (va, 4);       // shift the bits of the 16 src elements to right
            t1 = vandq_u8 (loset, va);    // obtain higher 4-bit of the src elements
            sh = vget_high_u8 (t1);
            sl = vget_low_u8 (t1);
            rh = vtbl2_u8 (mth, sh);
            rl = vtbl2_u8 (mth, sl);
            r2 = vcombine_u8 (rl, rh);
            
            //r2 = vqtbl1q_u8 (mth, t1);   // obtain products of the higher 4-bit
            r  = veorq_u8 (r, r2);         // obtain final result of src * multiplier
            //r = _mm_xor_si128 (r, _mm_shuffle_epi8 (mth, t1));
            va = vld1q_u8 (dptr);
            r = veorq_u8 (r, va);
            vst1q_u8 (dptr, r);
        }
        dptr += 16;
        sptr += 16;

    }
    return;
#elif defined(ARM_NEON64)
    uint8_t *sptr, *dptr, *top;
    sptr = src;
    dptr = dst;
    top  = src + bytes;

    uint8_t *bh, *bl;
    uint8x16_t mth, mtl, loset;

    if (multiplier != 1) {
        /* half tables only needed for multiplier != 1 */
        bh = (uint8_t*) galois_half_mult_table_high;
        bh += (multiplier << 4);
        bl = (uint8_t*) galois_half_mult_table_low;
        bl += (multiplier << 4);
        // read split tables as 128-bit values
        mth = vld1q_u8(bh);
        mtl = vld1q_u8(bl);
		loset = vdupq_n_u8(0x0f);
    }

    uint8x16_t va, vb, r, t1, r2;
    while (sptr < top)
    {

        if (sptr + 16 > top) {
            /* remaining data doesn't fit into __m128i, do not use SSE */
            for (i=0; i<top-sptr; i++) {
                if (multiplier == 1)
                    *(dptr+i) ^= *(sptr+i);
                else
                    *(dptr+i) ^= galois_mult_table[((*(sptr+i))<<GF_POWER) | multiplier];
            }
            break;
        }


        va = vld1q_u8 (sptr);
        if (multiplier == 1) {
            /* just XOR */
            vb = vld1q_u8 (dptr);
            vb = veorq_u8(va, vb);
            vst1q_u8 (dptr, vb);
        } else {
            /* use half tables */
            t1 = vandq_u8 (loset, va);    // obtain lower 4-bit of the 16 src elements
            r  = vqtbl1q_u8 (mtl, t1);    // obtain products of the lower 4-bit 
            va = vshrq_n_u8 (va, 4);       // shift the bits of the 16 src elements to right
            t1 = vandq_u8 (loset, va);    // obtain higher 4-bit of the src elements
            r2 = vqtbl1q_u8 (mth, t1);   // obtain products of the higher 4-bit
            r  = veorq_u8 (r, r2);         // obtain final result of src * multiplier
            //r = _mm_xor_si128 (r, _mm_shuffle_epi8 (mth, t1));
            va = vld1q_u8 (dptr);
            r = veorq_u8 (r, va);
            vst1q_u8 (dptr, r);
        }
        dptr += 16;
        sptr += 16;

    }
    return;
#else
    if (multiplier == 1) {
        for (i=0; i<bytes; i++)
            dst[i] ^= src[i];
        return;
    }

    for (i = 0; i < bytes; i++)
        dst[i] ^= galois_mult_table[(src[i]<<GF_POWER) | multiplier];
    return;
#endif
}

/*
 * Muliply a region of elements with multiplier. When NEON is available, use it
 */
void galois_multiply_region(uint8_t *src, uint8_t multiplier, int bytes)
{
    if (multiplier == 0) {
        memset(src, 0, sizeof(uint8_t)*bytes);
        return;
    } else if (multiplier == 1) {
        return;
    }
#if defined(ARM_NEON32)
    uint8_t *sptr, *top;
    sptr = src;
    top  = src + bytes;

    uint8_t *bh, *bl;
    uint8x8x2_t mth, mtl;
    uint8x16_t loset;
    /* half tables only needed for multiplier != 1 */
    bh = (uint8_t*) galois_half_mult_table_high;
    bh += (multiplier << 4);
    bl = (uint8_t*) galois_half_mult_table_low;
    bl += (multiplier << 4);
    // read split tables as 128-bit values
    mth.val[0] = vld1_u8(bh);
    mth.val[1] = vld1_u8(bh+8);
    mtl.val[0] = vld1_u8(bl);
    mtl.val[1] = vld1_u8(bl+8);
    loset = vdupq_n_u8(0x0f);

    uint8x16_t va, r, t1, r2;
    uint8x8_t sh, sl, rh, rl;
    while (sptr < top)
    {

        if (sptr + 16 > top) {
            /* remaining data doesn't fit into __m128i, do not use SSE */
            for (int i=0; i<top-sptr; i++)
                *(sptr+i) = galois_mult_table[((*(sptr+i))<<GF_POWER) | multiplier];
            break;
        }
        va = vld1q_u8 (sptr);
        t1 = vandq_u8 (loset, va);
        sh = vget_high_u8 (t1);
        sl = vget_low_u8 (t1);
        rh = vtbl2_u8 (mtl, sh);
        rl = vtbl2_u8 (mtl, sl);
        r  = vcombine_u8 (rl, rh);
        
        
        //r = vqtbl1q_u8 (mtl, t1);
        va = vshrq_n_u8 (va, 4);
        t1 = vandq_u8 (loset, va);
        sh = vget_high_u8 (t1);
        sl = vget_low_u8 (t1);
        rh = vtbl2_u8 (mth, sh);
        rl = vtbl2_u8 (mth, sl);
        r2 = vcombine_u8 (rl, rh);
        r  = veorq_u8 (r, r2);
        
        //r = veorq_u8 (r, vqtbl1q_u8 (mth, t1));
        vst1q_u8 (sptr, r);
        sptr += 16;

    }
    return;
#elif defined(ARM_NEON64)
    uint8_t *sptr, *top;
    sptr = src;
    top  = src + bytes;

    uint8_t *bh, *bl;
    uint8x16_t mth, mtl, loset;
    /* half tables only needed for multiplier != 1 */
    bh = (uint8_t*) galois_half_mult_table_high;
    bh += (multiplier << 4);
    bl = (uint8_t*) galois_half_mult_table_low;
    bl += (multiplier << 4);
    // read split tables as 128-bit values
    mth = vld1q_u8(bh);
    mtl = vld1q_u8(bl);
    loset = vdupq_n_u8(0x0f);

    uint8x16_t va, r, t1;
    while (sptr < top)
    {

        if (sptr + 16 > top) {
            /* remaining data doesn't fit into __m128i, do not use SSE */
            for (int i=0; i<top-sptr; i++)
                *(sptr+i) = galois_mult_table[((*(sptr+i))<<GF_POWER) | multiplier];
            break;
        }
        va = vld1q_u8 (sptr);
        t1 = vandq_u8 (loset, va);
        r = vqtbl1q_u8 (mtl, t1);
        va = vshrq_n_u8 (va, 4);
        t1 = vandq_u8 (loset, va);
        r = veorq_u8 (r, vqtbl1q_u8 (mth, t1));
        vst1q_u8 (sptr, r);
        sptr += 16;

    }
    return;
#else
    for (int i=0; i<bytes; i++)
        src[i] = galois_mult_table[((src[i])<<GF_POWER) | multiplier];
    return;
#endif
}
