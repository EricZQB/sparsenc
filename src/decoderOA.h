#ifndef OA_DECODER_H
#define OA_DECODER_H
#include "sparsenc.h"
struct running_matrix;
/*
 * OA (overlap aware) DECODING CONTEXT is used at each destination node
 * to record the status of the decoder. It is the core
 * information needed by the decoder
 */
struct decoding_context_OA
{
    // GNC context
    struct snc_context *sc;

    int aoh;                            // Allowed number of overhead packets (OHS)
    int finished;                       // an indicator tracking the finish of decoding
    int OA_ready;                       // the decoder has reached the necessary condition for OA decoding
    int local_DoF;                      // total DOF that have been received within generations, total_DoF==NUM_SRC, then decodable
    int global_DoF;                     // total true DoF that the receiver has received

    // Local decoding matrices
    struct running_matrix **Matrices;   //[CLASS_NUM] record running matrices of each class
    // Global decoding matrix
    GF_ELEMENT **JMBcoefficient;        //[NUM_SRC+OHS+CHECKS][NUM_PP];
    GF_ELEMENT **JMBmessage;            //[NUM_SRC+OHS+CHECKS][EXT_N];

    // the following two mappings are to record pivoting processings
    int *otoc_mapping;                  //[NUM_PP] record the mapping from original packet id to current column index
    int *ctoo_mapping;                  //[NUM_PP] record the mapping from current column index to the original packet id
    int inactives;                      // total number of inactivated packets among overlapping packets

    int overhead;                       // record how many packets have been received
    long long operations;               // record the number of computations used
};

struct decoding_context_OA *create_dec_context_OA(struct snc_parameter *sp, int aoh);
void process_packet_OA(struct decoding_context_OA *dec_ctx, struct snc_packet *pkt);
void free_dec_context_OA(struct decoding_context_OA *dec_ctx);

/**
 * File format to store ongoing decoding context
 *
 * snc_parameter
 * decoder_type
 * decoding_context_OA (excluding snc_context)
 *
 */
long save_dec_context_OA(struct decoding_context_OA *dec_ctx, const char *filepath);
struct decoding_context_OA *restore_dec_context_OA(const char *filepath);
#endif
