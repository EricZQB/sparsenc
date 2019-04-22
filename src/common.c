/*
 * Common utility functions used by many routines in the library.
 */
#include <stdint.h>
#include "common.h"
#include "galois.h"
int BALLOC = 500;

static int loglevel = 0;    // log level for the library
static int compare_int(const void *elem1, const void *elem2);
void set_loglevel(const char *level)
{
    if (strcmp(level, "TRACE") == 0)
        loglevel = TRACE;
    return;
}

int get_loglevel()
{
    return loglevel;
}

// check if an item is existed in an int array
int has_item(int array[], int item, int length)
{
    int index = -1;
    for (int i=0; i<length; i++) {
        if (item == array[i]) {
            index = i;
            break;
        }
    }
    return index;
}

void append_to_list(struct node_list *list, struct node *nd)
{
    if (list->first == NULL)
        list->first = list->last = nd;
    else {
        list->last->next = nd;
        list->last = nd;
    }
}

// Remove the first node whose data is equal to "data"
// Note: this function should only be used in applications
//       where nodes in the list have unique data
int remove_from_list(struct node_list *list, int data)
{
    struct node *prev = NULL;
    struct node *curr = list->first;
    while (curr != NULL) {
        if (curr->data == data) {
            // shorten list
            if (curr == list->first && curr == list->last) {            // list contains only one node
                list->first = list->last = NULL;
            } else if (curr == list->first && curr != list->last) {     // head node is to be removed
                list->first = curr->next;
            } else if (curr != list->first && curr == list->last) {     // tail node is to be removed
                list->last = prev;
                list->last->next = NULL;
            } else {
                prev->next = curr->next;
            }

            free(curr);
            curr = NULL;
            return 0;
        }
        prev = curr;
        curr = curr->next;
    }
    return -1;
}

int exist_in_list(struct node_list *list, int data)
{
    struct node *p = list->first;
    while (p != NULL) {
        if (p->data == data)
            return 1;
        p = p->next;
    }
    return 0;
}

// clear nodes in a list, but keep the list structure alive
void clear_list(struct node_list *list)
{
    struct node *nd = list->first;
    struct node *ndNext;
    while (nd != NULL) {
        ndNext = nd->next;
        free(nd);
        nd = ndNext;
    }
    list->first = list->last = NULL;
}

// Free a list, which include clear nodes in a list and free
// the list structure in the end.
void free_list(struct node_list *list)
{
    if (list != NULL) {
        clear_list(list);
        free(list);
    }
}

/**
 * Pack the small integers of equal length than len-bits in compact form in an array of bytes.
 * e.g. suppose all co's are <8, then every 3 bits represent an integer in the memory
 * pointed by coes
 * length - number of bits every co may occupy (Note: len <= 8)
 * i - the itndex of the such integer co in the memory from the beginning of coes[0]
 *     i=0,1,...
 *
 * Pay attention to the boundary of the byte array
 */

inline void pack_bits_in_byte_array(unsigned char *coes, int nbytes, unsigned char co, int len, int i)
{
    // Memory check
    if ((i - 1) * len >= 8 * nbytes) {
        printf("Error: accessing unauthorized memeory area.\n");
        exit(1);
    }
    int h = (len * i) % 8;          // position of the first bit of co in a byte
    int n = (len * i) / 8;          // index of byte where the first bit of co falls in
    // E.g.
    // coes:    [ ][ ][ ][ ][ ][ ][ ][ ]   [ ][ ][ ][ ][ ][ ][ ][ ]
    // co                          x  x     x  x  x  x
    // len=6, i=1: h=6, n=0
    //
    // coes:    [ ][ ][ ][ ][ ][ ][ ][ ]   [ ][ ][ ][ ][ ][ ][ ][ ]
    // co                                      x  x  x
    // len=3, i=3: h=1, n=1
    //

    unsigned char x, y, mask;
    // determine whether the integer occupies bits across bytes
    if (h + len <= 8) {
        // the i-th integer is in one byte
        x = co << (8 - (h + len)); //
        // clear the corresponding bits
        mask = ~((0xffU >> (8 - len)) << (8 - (h + len)));
        coes[n] &= mask;
        // set the corresponding bits
        coes[n] |= x;
    } else {
        // the i-th integer is across two bytes
        int lo = (h + len) % 8;     // number of left-over bits on the next byte
        // part 1: clear and set the corresponding bits on the tail of the first byte;
        x = co >> lo;                       // higher bits of co
        mask = ~(0xffU >> (8 - (len-lo)));
        coes[n] &= mask;
        coes[n] |= x;
        // part 2: clear and set the corresponding bits at the head of the next byte
        // if not reaching the boundary of the byte array.
        if (len * i > 8 * nbytes) {
            printf("pack_bits_in_byte_array(): reaching boundary.\n");
        } else {
            y = co << (8 - lo);             // lower bits of co (now at the head of a byte)
            mask = ~(0xffU << (8 - lo));
            coes[n+1] &= mask;              // clear the previous bits of coes[n+1]
            coes[n+1] |= y;
        }
    }
    return;
}

inline unsigned char read_bits_from_byte_array(unsigned char *coes, int nbytes, int len, int i)
{
    // Memory check
    if ((i - 1) * len >= 8 * nbytes) {
        printf("Error: accessing unauthorized memeory area.\n");
        exit(1);
    }
    int h = (len * i) % 8;          // position of the first bit of co in a byte
    int n = (len * i) / 8;          // index of byte where the first bit of co falls in
    unsigned char x, y, mask;
    unsigned char co = 0;
    if (h + len <=8) {
        mask = 0xffU >> (8 - len);         // len-1's mask
        mask = mask << (8 - (h + len));
        co = (mask & coes[n]) >> (8 - (h + len));
    } else {
        int lo = (h + len) % 8;         // number of bits on the "next" byte, so (len-lo) on the former
        x = coes[n] & (0xffU >> (8 - (len - lo)));  // bits at the tail of the first byte
        if (len * i <= 8 * nbytes) {
            y = coes[n+1] >> (8 - lo);              // bits at the head of the next byte
            co = (x << lo) | y;                     // combine the two parts
        } else {
            // Reaching the boundary, so only read bits in x
            printf("read_bits_from_byte_array(): reaching boundary.\n");
            co = (x << lo);                         // only read the higher bits of the intented len bits
        }
    }
    return co;
}

// A wrapper function of multiply_add_region for GF(2^2), ..., GF(2^7)
// nelem - number of GF(4), GF(8), GF(16), ..., GF(128) elements
// nbytes - number of bytes pointed by dst and src
void galois2n_multiply_add_region(GF_ELEMENT *dst, GF_ELEMENT *src, GF_ELEMENT multiplier, int gfpower, int nelem, int nbytes)
{
    if (multiplier == 0) {
        // add nothing
        return;
    }
    int i, j;
    if (multiplier == 1) {
        // add is XOR for all GF(2) extension fields
        for (i=0; i<nbytes; i++)
            dst[i] ^= src[i];
        return;
    }
    // Complication is here
    for (j=0; j<nelem; j++) {
        GF_ELEMENT dst_elem = read_bits_from_byte_array(dst, nbytes, gfpower, j);
        GF_ELEMENT new_elem = galois_add(dst_elem, galois_multiply(read_bits_from_byte_array(src, nbytes, gfpower, j), multiplier));
        pack_bits_in_byte_array(dst, nbytes, new_elem, gfpower, j);
    }
}

/**
 * Get/set the i-th bit from a sequence of bytes pointed
 * by coes. The indices of bits are as following:
 *
 *    [7|6|5|4|3|2|1|0]   [15|14|13|12|11|10|9|8]   ...
 *
 * It's caller's responsibility to ensure that ceil(max(i)/8) elements
 * are allocated in the memory pointed by coes.
 */
inline unsigned char get_bit_in_array(unsigned char *coes, int i)
{
    unsigned char co = coes[i/8];
    unsigned char mask = 0x1 << (i % 8);
    return ((mask & co) == mask);
}

inline void set_bit_in_array(unsigned char *coes, int i)
{
    coes[i/8] |= (0x1 << (i % 8));
    return;
}

// Build subgeneration neighbors list for each packet according to subgeneration grouping scheme
ID_list **build_subgen_nbr_list(struct snc_context *sc)
{
    int numpp = sc->snum + sc->cnum;
    ID_list **gene_nbr = malloc(sizeof(ID_list*) * numpp);
    int i, j;
    for (i=0; i<numpp; i++) {
        gene_nbr[i] = calloc(1, sizeof(ID_list));
    }
    for (i=0; i<sc->gnum; i++) {
        for (j=0; j<sc->params.size_g; j++) {
            int pktid = sc->gene[i]->pktid[j];
            if (!exist_in_list(gene_nbr[pktid], i)) {
                ID *nb = calloc(1, sizeof(ID));
                if (nb == NULL)
                    return NULL;
                nb->data = i;
                nb->ce   = 1;
                nb->next = NULL;
                append_to_list(gene_nbr[pktid], nb);
            }
        }
    }
    return gene_nbr;
}

void free_subgen_nbr_list(struct snc_context *sc, ID_list **gene_nbr)
{
    if (gene_nbr == NULL)
        return;
    int numpp = sc->snum + sc->cnum;
    for (int i=0; i<numpp; i++) {
        if (gene_nbr[i] != NULL) {
            free_list(gene_nbr[i]);
            gene_nbr[i] = NULL;
        }
    }
    return;
}

// generate a number of n<ub unique random numbers within the range of [0, ub-1]
// using Fisher-Yates shuffle method
void get_random_unique_numbers(int ids[], int n, int ub)
{
	int init_array[ub];
	int i, j;
	for (i=0; i<ub; i++)
		init_array[i] = i;

	// randomly shuffle the init_array
	for (i=ub-1; i>=1; i--) {
		int rd = genrand_int32() % (i+1);
		//int rd = gsl_rng_uniform_int(r, i+1);
		int tmp = init_array[rd];
		init_array[rd] = init_array[i];
		init_array[i] = tmp;
	}

    // sort the obtained unique random numbers so that coding coefficients corresponding
    // to packets are stored in the ascending order (to simplify decoder implementation)
    qsort(init_array, n, sizeof(int), compare_int);
    memcpy(ids, init_array, n*sizeof(int));
	//for (j=0; j<n; j++)
	//	ids[j] = init_array[j];
}

static int compare_int(const void *elem1, const void *elem2)
{
    int a = * ((int *) elem1);
    int b = * ((int *) elem2);
    if (a < b)
        return -1;
    if (a > b)
        return 1;
    return 0;
}
/*
 * Swap two continuous memory blocks
 */
/*
   void swap_memory(uint8_t *a1, uint8_t *a2, int bytes)
   {
   int i;
#if defined(INTEL_SSSE3)
uint8_t *sptr, *dptr, *top;
sptr = a1;
dptr = a2;
top  = a1 + bytes;

__m128i va, vb, r, t1;
while (sptr < top)
{
if (sptr + 16 > top) {
// remaining data doesn't fit into __m128i, do not use SSE
for (i=0; i<top-sptr; i++) {
uint8_t temp = *(dptr+i);
 *(dptr+i) = *(sptr+i);
 *(sptr+i) = temp;
 }
 break;
 }
 va = _mm_loadu_si128 ((__m128i *)(sptr));
 vb = _mm_loadu_si128 ((__m128i *)(dptr));
 _mm_storeu_si128 ((__m128i *)(dptr), va);
 _mm_storeu_si128 ((__m128i *)(sptr), vb);
 dptr += 16;
 sptr += 16;
 }
 return;
#else
for (i = 0; i < bytes; i++) {
uint8_t temp = a2[i];
a2[i] = a1[i];
a1[i] = temp;
}
return;
#endif
}
*/

/*
static unsigned long int next = 1;
int snc_rand(void)
{
    next = next * 1103515245 + 12345;
    return (unsigned int)(next/65536) % 32768;
}

void snc_srand(unsigned int seed)
{
    next = seed;
}
*/
