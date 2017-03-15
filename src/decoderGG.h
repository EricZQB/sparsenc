#ifndef NGG_DECODER_H
#define NGG_DECODER_H
#include "sparsenc.h"

#define FB_THOLD    1
typedef struct node      ID;
typedef struct node_list ID_list;

// Store in-processing matrices of the decoder
struct running_matrix;

struct decoding_context_GG {
    struct snc_context  *sc;            // The file information
    /********************************
     * Used in decoding LDPC pre-code
     ********************************/
    GF_ELEMENT **evolving_checks;       // Evolving check packets during precode decoding
    // We cannot change decoded packets in file_context directly during iterative decoding,
    // because those packets may be used in decoding other generations later
    int     *check_degrees;             // The number of unknown(undeocded) source neighbors of each check packet

    /******************************
     * Used in decoding generations
     ******************************/
    int finished;
    int decoded;                        // record how many packets have been decoded
    int originals;                      // record how many source packets are decoded
    struct running_matrix **Matrices;   // record running matrices of each class
    ID_list *recent;                    // record most recently decoded packets IDs
    /*******************************************
     * Used if feedback to encoder is allowed
     ******************************************/
    int grecent[FB_THOLD];              // gids of recent FB_THOLD decoded generations (it would wrap around)
    int newgpos;                        // pos pointer in grecent where newly decoded gid is stored
    int grcount;                        // the number of decoded generations counted from last feedback
    /*******************************************
     * Record for anaylzing decoding performance
     * Not actually needed in implementation
     *******************************************/
    long long operations;               // record the number of computations used
    int overhead;                       // record how many packets have been received
    long long ops1, ops2;               // splitted decoding operations
                                        // ops1 - operations invoked by Gaussian elimination
                                        // ops2 - operations invoked by substitutions
};

/**
 * Create the GG decode context
 * Return Values
 *   On success - return 0
 *   Otherwise  - return -1
 */
struct decoding_context_GG *create_dec_context_GG(struct snc_parameters *sp);
void free_dec_context_GG(struct decoding_context_GG *dec_ctx);
void process_packet_GG(struct decoding_context_GG *dec_ctx, struct snc_packet *pkt);

/**
 * File format to store ongoing decoding context
 *
 * snc_parameter
 * decoder_type
 * decoding_context_GG (excluding snc_context)
 *
 */
long save_dec_context_GG(struct decoding_context_GG *dec_ctx, const char *filepath);
struct decoding_context_GG *restore_dec_context_GG(const char *filepath);
#endif /* GG_DECODER */
