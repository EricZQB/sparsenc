/*
 * Common utility functions used by many routines in the library.
 */
#include <stdint.h>
#include "common.h"
static int loglevel = 0;    // log level for the library
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
