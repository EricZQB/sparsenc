/*--------------------- common.h ------------------------
 *  Internal header file
 *------------------------------------------------------*/
#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "bipartite.h"
#include "sparsenc.h"

/* log levels */
#define TRACE       5

#define ALIGN(a, b) ((a) % (b) == 0 ? (a)/(b) : (a)/(b) + 1)
#define RESIDUAL(a, b) ((b) * ALIGN((a), (b)) - (a))

#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
typedef struct node      ID;
typedef struct node_list ID_list;

/**
 * Source packets are grouped into subsets, referred to as
 * subgenerations in this library to emphasize the 'sub-' concept
 * in sparse code. (Subsets are also referred to as batches, classes,
 * chunks, segments, or generations in many other network coding literature).
 **/
struct subgeneration {
    int gid;
    int *pktid;                 // SIZE_G source packet IDs
};

/**
 * Definition of snc_context
 **/
struct snc_context {
    struct  snc_parameters    params;
    int                       snum;     // Number of source packets splitted
    int                       cnum;     // Number of parity-checks(cnum ~= snum * pcrate)
    int                       gnum;     // Number of subgenerations
    struct  subgeneration   **gene;     // array of pointers each points to a subgeneration.
    struct  bipartite_graph  *graph;
    GF_ELEMENT              **pp;       // Pointers to precoded source packets
    int                      *nccount;  // Count of coded packets generated from each subgeneration
    int                       count;    // Count of total coded packets generated
};


/*
 * Buffer for storing SNC packets (for recoding)
 *
 * Buffer size specifies how many packets are saved for
 * each subgeneration. "FIFO" strategy is used when buffer
 * size is reached; the oldest buffered packet will be
 * discarded when a subgeneration buffer is full while a new
 * packet belonging to the subgeneration arrives.
 *
 * Buffer data structure
 *
 * gbuf --> gbuf[0]
 *                     snc_packet  snc_packet ...
 *          gbuf[1]          ^            ^
 *                           |            |
 *          gbuf[2] --> gbuf[2][0]   gbuf[2][1] ....
 *            .
 *            .
 *            .
 *
 * Systematic packet buffer (systematic code)
 *        snc_packet
 *             ^
 *             |
 * sysbuf[0] sysbuf[1] sysbuf[2] ... 
 */
struct snc_buffer {
    struct snc_parameters  params;      // Meta info of the code
    int                    snum;        // Number of source packets
    int                    cnum;        // Number of parity-check packets
    int                    gnum;        // Number of subgenerations
    int                    size;        // Number of bufferred packets of each subgeneration
    int                    nemp;        // Number of non-empty subgeneration buffers
    struct snc_packet   ***gbuf;        // Pointers to subgeneration buffers
    int                   *nc;          // Number of currently buffered packets of each generation
    int                   *pn;          // Positions to store next packet of each subgeneration
    int                   *nsched;      // Number of scheduled times of each subgeneration
    // Use if code in buffer is systematic
    struct snc_packet    **sysbuf;      // Buffered uncoded packet (needed for systematic code)
    int                    spn;         // Position to store next systematic packet in sysbuf
    int                    sysnum;      // number of buffered systematic packet
    int                    newsys;      // A new systematic code has been received, schedule it
                                        // -1 if not received, >=0 indicates its position in the sysbuf
    int                    sysptr;      // pointer of already scheduled systematic packet
    int                    sys_sched;   // scheduled systematic packet index in sysbuf
};


// BATS buffer of fixed size
// Received packets from upstream are buffered in a first-in-first-out manner. 
// If the buffer is full, the oldest buffered packet would be discarded to store 
// a new received packet. We assume that the batch size is much larger than the 
// buffer size. Therefore, the buffered packets may belong to at most two different 
// batches, which happens when the buffer starts to receive the first several packets 
// belonging to a new batch. When there are two different batches in the buffer, 
// recoded packets are still generated from the older batch until all of its packets 
// are discarded to accommodate received packets of the new batch. The current batch 
// from which recoded packets are generated is referred to as the \textit{sending} batch 
// of the buffer, and the \textit{receiving} batch is the batch the latest received 
// packet belongs to. Clearly, the sending and receiving batches are the same if there 
// are only one batch in the buffer.
struct snc_buffer_bats {
    struct snc_parameters   params;             // pointer to the parameter of the BATS code 
    struct snc_packet     **srbuf;
    int                     bufsize;            // size of buffer
    int                     sbatchid;           // current sending batch
    int                     s_first;            // start pos index of sending buffer
    int                     r_last;             // end pos index of receiving buffer
                                                // if ((r_last+1) % bufsize == s_start), discard old pkt
};

/* Row vector of a matrix */
struct row_vector
{
    int len;            // length of the row
    GF_ELEMENT *elem;   // elements of the row
};


/* common.c */
void set_loglevel(const char *level);
int get_loglevel();
int has_item(int array[], int item, int length);
void append_to_list(struct node_list *list, struct node *nd);
int remove_from_list(struct node_list *list, int data);
int exist_in_list(struct node_list *list, int data);
void clear_list(struct node_list *list);
void free_list(struct node_list *list);
unsigned char get_bit_in_array(unsigned char *coes, int i);
void set_bit_in_array(unsigned char *coes, int i);
ID_list **build_subgen_nbr_list(struct snc_context *sc);
void free_subgen_nbr_list(struct snc_context *sc, ID_list **gene_nbr);
void get_random_unique_numbers(int ids[], int n, int ub);
//int snc_rand(void);
//void snc_srand(unsigned int seed);
// mt19937ar.c
void init_genrand(unsigned long s);
unsigned long genrand_int32(void);
#endif /* COMMON_H */
