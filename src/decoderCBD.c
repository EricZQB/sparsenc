/********************************************************************
 *                       Compact Band Decoder
 *
 * The decoder only applies to band code.
 *
 * Unlike regular band decoder (decoderBD.c), CBD decoder is
 * compact in coefficient matrix storage. Only coefficients in the band
 * (which are nonzeros) are stored. A price to pay is that pivoting
 * cannot be performed due to the limited random access and row/col
 * manipulation capability of using compact row vectors.
 ********************************************************************/
#include "common.h"
#include "galois.h"
#include "decoderCBD.h"
static int process_vector_CBD(struct decoding_context_CBD *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message);
static int apply_parity_check_matrix(struct decoding_context_CBD *dec_ctx);
static void finish_recovering_CBD(struct decoding_context_CBD *dec_ctx);

// create decoding context for band decoder
struct decoding_context_CBD *create_dec_context_CBD(struct snc_parameters *sp)
{
    static char fname[] = "snc_create_dec_context_CBD";
    int i, j, k;

    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    int niv = 0;
    if (sp->type != BAND_SNC) {
        fprintf(stdout, "WARNING: Band decoder only applies to band GNC code. Fallback to naive mode for non-band codes.\n");
        niv = 1;
    }

    struct decoding_context_CBD *dec_ctx = malloc(sizeof(struct decoding_context_CBD));
    if (dec_ctx == NULL) {
        fprintf(stderr, "%s: malloc decoding context CBD failed\n", fname);
        return NULL;
    }

    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }

    dec_ctx->sc = sc;

    dec_ctx->finished     = 0;
    dec_ctx->DoF          = 0;
    dec_ctx->de_precode   = 0;
    dec_ctx->naive        = niv;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    dec_ctx->row = (struct row_vector **) calloc(numpp, sizeof(struct row_vector *));
    if (dec_ctx->row == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->row failed\n", fname);
        goto AllocError;
    }
    dec_ctx->message = calloc(numpp, sizeof(GF_ELEMENT*));
    if (dec_ctx->message == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->message failed\n", fname);
        goto AllocError;
    }
    for (i=0; i<numpp; i++) {
        dec_ctx->message[i] = calloc(pktsize, sizeof(GF_ELEMENT));
        if (dec_ctx->message[i] == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->message[%d] failed\n", fname, i);
            goto AllocError;
        }
    }

    dec_ctx->overhead     = 0;
    dec_ctx->operations   = 0;
    dec_ctx->ops1 = 0;             // operations of forward sub
    dec_ctx->ops2 = 0;             // operations of applying precode
    dec_ctx->ops3 = 0;             // operations of backward sub
    return dec_ctx;

AllocError:
    free_dec_context_CBD(dec_ctx);
    return NULL;
}

/*
 * Note: throughout the packet collecting process, the decoding matrix
 * is maintained an upper triangular form.
 */
void process_packet_CBD(struct decoding_context_CBD *dec_ctx, struct snc_packet *pkt)
{
    static char fname[] = "snc_process_packet_CBD";
    dec_ctx->overhead += 1;
    int i, j, k;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // transform GNC encoding vector to full length
    GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));
    if (ces == NULL)
        fprintf(stderr, "%s: calloc ces failed\n", fname);

    // Pay attention to systematic packets. It is easy to handle in CBD decoder. Let's 
    // just transform it to full-length singleton vector.
    if (pkt->gid == -1 && pkt->ucid == -1) {
        fprintf(stderr, "%s: pkt's gid is -1 but ucid is not valid\n", fname);
    } else if (pkt->gid == -1 && pkt->ucid >= 0) {
        ces[pkt->ucid] = 1;
    } else {
        // This is normal GNC packet
        for (i=0; i<gensize; i++) {
            int index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
            if (dec_ctx->sc->params.bnc) {
                ces[index] = get_bit_in_array(pkt->coes, i);
            } else {
                ces[index] = pkt->coes[i];
            }
        }
    }

    /* Process full-length encoding vector against decoding matrix */
    int lastDoF = dec_ctx->DoF;
    int pivot = process_vector_CBD(dec_ctx, ces, pkt->syms);
    free(ces);
    ces = NULL;
    if (get_loglevel() == TRACE) 
        printf("received %d DoF: %d\n", dec_ctx->overhead, dec_ctx->DoF-lastDoF);
    // If the number of received DoF is equal to NUM_SRC, apply the parity-check matrix.
    // The messages corresponding to rows of parity-check matrix are all-zero.
    if (dec_ctx->DoF == dec_ctx->sc->snum) {
        dec_ctx->de_precode = 1;    /*Mark de_precode before applying precode matrix*/
        int missing_DoF = apply_parity_check_matrix(dec_ctx);
        if (get_loglevel() == TRACE)
            printf("After applying the parity-check matrix, %d DoF are missing.\n", missing_DoF);
        dec_ctx->DoF = numpp - missing_DoF;
    }

    if (dec_ctx->DoF == dec_ctx->sc->snum + dec_ctx->sc->cnum) {
        finish_recovering_CBD(dec_ctx);
    }
}

/* Process a full row vector against CBD decoding matrix */
static int process_vector_CBD(struct decoding_context_CBD *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message)
{
    static char fname[] = "process_vector_CBD";
    int i, j, k;
    int pivot = -1;
    int pivotfound = 0;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    int rowop = 0;
    for (i=0; i<numpp; i++) {
        if (vector[i] != 0) {
            if (dec_ctx->row[i] != NULL) {
                /* There is a valid row saved for pivot-i, process against it */
                assert(dec_ctx->row[i]->elem[0]);
                quotient = galois_divide(vector[i], dec_ctx->row[i]->elem[0]);
                galois_multiply_add_region(&(vector[i]), dec_ctx->row[i]->elem, quotient, dec_ctx->row[i]->len);
                galois_multiply_add_region(message, dec_ctx->message[i], quotient, pktsize);
                dec_ctx->operations += 1 + dec_ctx->row[i]->len + pktsize;
                if (!dec_ctx->de_precode) {
                    dec_ctx->ops1 += 1 + dec_ctx->row[i]->len + pktsize;
                } else {
                    dec_ctx->ops2 += 1 + dec_ctx->row[i]->len + pktsize;
                }
                rowop += 1;
            } else {
                pivotfound = 1;
                pivot = i;
                break;
            }
        }
    }

    if (pivotfound == 1) {
        /* Save it to the corresponding row */
        dec_ctx->row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
        if (dec_ctx->row[pivot] == NULL)
            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
        int len;
        if (!dec_ctx->de_precode && !dec_ctx->naive) {
            /* before de_precode every row is no more than gensize-width */
            len = numpp - pivot > gensize ? gensize : numpp - pivot;
        } else {
            /* row bandwidth is indetermined, so being conservative here */
            len = numpp - pivot;
        }
        dec_ctx->row[pivot]->len = len;
        dec_ctx->row[pivot]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
        if (dec_ctx->row[pivot]->elem == NULL)
            fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, pivot);
        memcpy(dec_ctx->row[pivot]->elem, &(vector[pivot]), len*sizeof(GF_ELEMENT));
        assert(dec_ctx->row[pivot]->elem[0]);
        memcpy(dec_ctx->message[pivot], message,  pktsize*sizeof(GF_ELEMENT));
        if (get_loglevel() == TRACE) 
            printf("received-DoF %d new-DoF %d row_ops: %d\n", dec_ctx->DoF, pivot, rowop);
        dec_ctx->DoF += 1;
    }
    return pivot;
}

// Apply the parity-check matrix to the decoding matrix
static int apply_parity_check_matrix(struct decoding_context_CBD *dec_ctx)
{
    static char fname[] = "apply_parity_check_matrix";
    int i, j, k;
    int num_of_new_DoF = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
    GF_ELEMENT *ces = malloc(numpp*sizeof(GF_ELEMENT));
    GF_ELEMENT *msg = malloc(pktsize*sizeof(GF_ELEMENT));
    int p = 0;          // index pointer to the parity-check vector that is to be copyed
    for (int p=0; p<dec_ctx->sc->cnum; p++) {
        memset(ces, 0, numpp*sizeof(GF_ELEMENT));
        memset(msg, 0, pktsize*sizeof(GF_ELEMENT));
        /* Set the coding vector according to parity-check bits */
        NBR_node *varnode = dec_ctx->sc->graph->l_nbrs_of_r[p]->first;
        while (varnode != NULL) {
            ces[varnode->data] = varnode->ce;
            varnode = varnode->next;
        }
        ces[dec_ctx->sc->snum+p] = 1;
        int pivot = process_vector_CBD(dec_ctx, ces, msg);
    }
    free(ces);
    free(msg);

    /* Count available innovative rows */
    int missing_DoF = 0;
    for (i=0; i<numpp; i++) {
        if (dec_ctx->row[i] == NULL)
            missing_DoF++;
        else if (dec_ctx->row[i]->elem[0] ==0 && get_loglevel() == TRACE) {
            printf("%s: row[%d]->elem[0] is 0\n", fname, i);
        }
    }
    return missing_DoF;
}

/**
 * Finish CBD decoding
 * This routine converts decoding matrix from upper triangular
 * form to diagonal.
 */
static void finish_recovering_CBD(struct decoding_context_CBD *dec_ctx)
{
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    int i, j;
    int len;
    GF_ELEMENT quotient;
    for (i=numpp-1; i>=0; i--) {
        /* eliminate all nonzeros above diagonal elements from right to left*/
        for (j=0; j<i; j++) {
            len = dec_ctx->row[j]->len;
            if (j+len <= i || dec_ctx->row[j]->elem[i-j] == 0)
                continue;
            assert(dec_ctx->row[i]->elem[0]);
            quotient = galois_divide(dec_ctx->row[j]->elem[i-j], dec_ctx->row[i]->elem[0]);
            galois_multiply_add_region(dec_ctx->message[j], dec_ctx->message[i], quotient, pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->ops3 += (pktsize + 1);
            dec_ctx->row[j]->elem[i-j] = 0;
        }
        /* convert diagonal to 1*/
        if (dec_ctx->row[i]->elem[0] != 1) {
            galois_multiply_region(dec_ctx->message[i], galois_divide(1, dec_ctx->row[i]->elem[0]), pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->ops3 += (pktsize + 1);
            dec_ctx->row[i]->elem[0] = 1;
        }
        /* save decoded packet */
        dec_ctx->sc->pp[i] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->sc->pp[i], dec_ctx->message[i], pktsize*sizeof(GF_ELEMENT));
    }
    dec_ctx->finished = 1;
    if (get_loglevel() == TRACE) {
        int snum = dec_ctx->sc->snum;
        printf("Splitted operations: %f %f %f\n", (double) dec_ctx->ops1/snum/pktsize,
                                                  (double) dec_ctx->ops2/snum/pktsize,
                                                  (double) dec_ctx->ops3/snum/pktsize);
    }
}

void free_dec_context_CBD(struct decoding_context_CBD *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    if (dec_ctx->sc != NULL)
        snc_free_enc_context(dec_ctx->sc);
    if (dec_ctx->row != NULL) {
        for (int i=dec_ctx->sc->snum+dec_ctx->sc->cnum-1; i>=0; i--) {
            if (dec_ctx->row[i] != NULL) {
                if (dec_ctx->row[i]->elem != NULL)
                    free(dec_ctx->row[i]->elem);
                free(dec_ctx->row[i]);
            }
        }
        free(dec_ctx->row);
    }
    if (dec_ctx->message != NULL) {
        for (int i=dec_ctx->sc->snum+dec_ctx->sc->cnum-1; i>=0; i--) {
            if (dec_ctx->message[i] != NULL)
                free(dec_ctx->message[i]);
        }
        free(dec_ctx->message);
    }
    free(dec_ctx);
    dec_ctx = NULL;
    return;
}

/**
 * Save a decoding context to a file
 * Return values:
 *   On success: bytes written
 *   On error: -1
 */
long save_dec_context_CBD(struct decoding_context_CBD *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = CBD_DECODER;
    FILE *fp;
    if ((fp = fopen(filepath, "w")) == NULL) {
        fprintf(stderr, "Cannot open %s to save decoding context\n", filepath);
        return (-1);
    }
    // Write snc params
    filesize += fwrite(&dec_ctx->sc->params, sizeof(struct snc_parameters), 1, fp);
    // Write decoder type
    filesize += fwrite(&(d_type), sizeof(int), 1, fp);
    // Write decoding context
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->DoF, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->de_precode, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->naive, sizeof(int), 1, fp);
    // Decoder matrix rows
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    for (int i=0; i<numpp; i++) {
        int rowlen = dec_ctx->row[i] == NULL ? 0 : dec_ctx->row[i]->len;
        filesize += fwrite(&rowlen, sizeof(int), 1, fp);
        if (rowlen != 0) {
            filesize += fwrite(dec_ctx->row[i]->elem, sizeof(GF_ELEMENT), rowlen, fp);
            filesize += fwrite(dec_ctx->message[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    /*performance index*/
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

/*
 * Restore a decoding context from a file
 */
struct decoding_context_CBD *restore_dec_context_CBD(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    // Create a fresh decoding context
    struct decoding_context_CBD *dec_ctx = create_dec_context_CBD(&sp);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_CBD failed\n");
        return NULL;
    }
    // Restore decoding context from file
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field

    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->DoF, sizeof(int), 1, fp);
    fread(&dec_ctx->de_precode, sizeof(int), 1, fp);
    fread(&dec_ctx->naive, sizeof(int), 1, fp);
    // Decoder matrix rows
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    for (int i=0; i<numpp; i++) {
        int rowlen = 0;
        fread(&rowlen, sizeof(int), 1, fp);
        if (rowlen != 0) {
            // There is an non-NULL row
            dec_ctx->row[i] = (struct row_vector*) malloc(sizeof(struct row_vector));
            if (dec_ctx->row[i] == NULL) {
                free_dec_context_CBD(dec_ctx);
                return NULL;
            }
            dec_ctx->row[i]->len = rowlen;
            dec_ctx->row[i]->elem = (GF_ELEMENT *) calloc(rowlen, sizeof(GF_ELEMENT));
            if (dec_ctx->row[i]->elem == NULL) {
                free_dec_context_CBD(dec_ctx);
                return NULL;
            }
            fread(dec_ctx->row[i]->elem, sizeof(GF_ELEMENT), rowlen, fp);
            fread(dec_ctx->message[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    /*performance index*/
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
