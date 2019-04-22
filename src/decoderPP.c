/********************************************************************
 *                       Perpetual Code Decoder
 *
 * The decoder only applies to perpetual code which has band structure
 * with wrap-around in some rows.
 ********************************************************************/
#include "common.h"
#include "galois.h"
#include "decoderPP.h"
static int process_vector_PP(struct decoding_context_PP *dec_ctx, GF_ELEMENT *vector, GF_ELEMENT *message);
static void finish_recovering_PP(struct decoding_context_PP *dec_ctx);

// create decoding context for perpetual decoder
struct decoding_context_PP *create_dec_context_PP(struct snc_parameters *sp)
{
    static char fname[] = "snc_create_dec_context_CBD";
    int i, j, k;

    if (sp->type != WINDWRAP_SNC || sp->size_c != 0) {
        fprintf(stdout, "WARNING: PP decoder only applies to perpetual codes.\n");
        exit(1);
    }

    struct decoding_context_PP *dec_ctx = malloc(sizeof(struct decoding_context_PP));
    if (dec_ctx == NULL) {
        fprintf(stderr, "%s: malloc decoding context PP failed\n", fname);
        return NULL;
    }

    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }

    dec_ctx->sc = sc;

    dec_ctx->stage        = FORWARD;
    dec_ctx->pivots       = 0;
    dec_ctx->finished     = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum;

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
    return dec_ctx;

AllocError:
    free_dec_context_PP(dec_ctx);
    return NULL;
}

void process_packet_PP(struct decoding_context_PP *dec_ctx, struct snc_packet *pkt)
{
    static char fname[] = "snc_process_packet_PP";
    dec_ctx->overhead += 1;

    
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum;

    // start processing
    int i, j, k;
    GF_ELEMENT quotient;
    if (dec_ctx->stage == FORWARD) {
        // transform GNC encoding vector to full length (gensize) in case it is GF(2) and therefore was compressed
        GF_ELEMENT *ces0 = calloc(gensize, sizeof(GF_ELEMENT));
        if (ces0 == NULL)
            fprintf(stderr, "%s: calloc ces0 failed\n", fname);
        if (dec_ctx->sc->params.gfpower ==1) {
            for (i=0; i<gensize; i++)
                ces0[i] = get_bit_in_array(pkt->coes, i);
        } else if (dec_ctx->sc->params.gfpower == 8) {
            memcpy(ces0, pkt->coes, gensize*sizeof(GF_ELEMENT));
        } else {
            for (i=0; i<gensize; i++)
                ces0[i] = read_bits_from_byte_array(pkt->coes, dec_ctx->sc->params.size_g, dec_ctx->sc->params.gfpower, i);
        }
        int pivot = dec_ctx->sc->gene[pkt->gid]->pktid[0];// By default, the coding coefficient of the pivot candidate is pkt->coes[0]
        int shift = 0;
        while (ces0[shift] == 0) {
            shift += 1;
            if (shift == gensize) {
                free(ces0);
                return;         // pkt->coes is all zero, so this is a useless packet, just return.
            }
        }
        pivot = (pivot + shift) % numpp;
        GF_ELEMENT *ces_tmp = calloc(gensize, sizeof(GF_ELEMENT));
        memcpy(ces_tmp, &(ces0[shift]), (gensize-shift)*sizeof(GF_ELEMENT));  // temporary place for multiply-add with existing rows
        // There is already a row with the same pivot in decoding matrix
        int rowlen = gensize - shift;
        while (dec_ctx->row[pivot] != NULL) {
            quotient = galois_divide(ces_tmp[0], dec_ctx->row[pivot]->elem[0]);
            galois_multiply_add_region(ces_tmp, dec_ctx->row[pivot]->elem, quotient, dec_ctx->row[pivot]->len);
            galois_multiply_add_region(pkt->syms, dec_ctx->message[pivot], quotient, pktsize);
            int newlen = rowlen > dec_ctx->row[pivot]->len ? rowlen : dec_ctx->row[pivot]->len;  // new length of the vector after processed
            rowlen = newlen - 1;   // the first element has been reduced to 0, so omit it
            memset(ces0, 0, sizeof(GF_ELEMENT)*gensize);
            memcpy(ces0, &(ces_tmp[1]), rowlen*sizeof(GF_ELEMENT));    // copy resultant vector back to ces0
            dec_ctx->operations += 1 + dec_ctx->row[pivot]->len + pktsize;
            shift = 0;
            while (ces0[shift] == 0) {
                shift += 1;
                if (shift == newlen) {
                    free(ces0);
                    free(ces_tmp);
                    return;         // pkt->coes is reduced to zero, so this is a useless packet, just return.
                }
            }
            pivot = (pivot + 1 + shift) % numpp;
            memset(ces_tmp, 0, sizeof(GF_ELEMENT)*gensize);
            memcpy(ces_tmp, &(ces0[shift]), (rowlen-shift)*sizeof(GF_ELEMENT));
            rowlen = rowlen - shift;
        }
        // Save the resultant row
        dec_ctx->row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
        if (dec_ctx->row[pivot] == NULL)
            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
        int len = rowlen;
        dec_ctx->row[pivot]->len = len;
        dec_ctx->row[pivot]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
        if (dec_ctx->row[pivot]->elem == NULL)
            fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, pivot);
        memcpy(dec_ctx->row[pivot]->elem, ces_tmp, len*sizeof(GF_ELEMENT));
        assert(dec_ctx->row[pivot]->elem[0]);
        memcpy(dec_ctx->message[pivot], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
        dec_ctx->pivots += 1;
        if (get_loglevel() == TRACE) 
            printf("pivot-candidates %d received %d\n", dec_ctx->pivots, dec_ctx->overhead);
        free(ces0);
        free(ces_tmp);
    } else if (dec_ctx->stage == FINALFORWARD) {
        // Previous final forward was not successful, and therefore it receives more packets to fill in the decoding matrix
        // Now always convert encoding vector to full length (numpp)
        GF_ELEMENT *ces1 = calloc(numpp, sizeof(GF_ELEMENT));
        if (ces1 == NULL)
            fprintf(stderr, "%s: calloc ces1 failed\n", fname);
        for (i=0; i<gensize; i++) {
            int index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
            if (dec_ctx->sc->params.gfpower == 1) {
                ces1[index] = get_bit_in_array(pkt->coes, i);
            } else if (dec_ctx->sc->params.gfpower == 8) {
                ces1[index] = pkt->coes[i];
            } else {
                ces1[index] = read_bits_from_byte_array(pkt->coes, dec_ctx->sc->params.size_g, dec_ctx->sc->params.gfpower, i);
            }
        }
        // Process the full length vector against existing rows
        for (k=0; k<numpp; k++) {
            if (ces1[k] != 0) {
                if (dec_ctx->row[k] != NULL) {
                    assert(dec_ctx->row[k]->elem[0]);
                    quotient = galois_divide(ces1[k], dec_ctx->row[k]->elem[0]);
                    galois_multiply_add_region(&(ces1[k]), dec_ctx->row[k]->elem, quotient, dec_ctx->row[k]->len);
                    galois_multiply_add_region(pkt->syms, dec_ctx->message[k], quotient, pktsize);
                    dec_ctx->operations += 1 + dec_ctx->row[k]->len + pktsize;
                } else {
                    // a valid pivot found, store it back to decoding matrix
                    dec_ctx->row[k] = (struct row_vector*) malloc(sizeof(struct row_vector));
                    if (dec_ctx->row[k] == NULL)
                        fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, k);
                    int len = numpp - k;
                    dec_ctx->row[k]->len = len;
                    dec_ctx->row[k]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
                    if (dec_ctx->row[k]->elem == NULL)
                        fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, k);
                    memcpy(dec_ctx->row[k]->elem, &(ces1[k]), len*sizeof(GF_ELEMENT));
                    assert(dec_ctx->row[k]->elem[0]);
                    memcpy(dec_ctx->message[k], pkt->syms, pktsize*sizeof(GF_ELEMENT));
                    dec_ctx->pivots += 1;
                    break;
                }
            }
        }
        free(ces1);
        if (dec_ctx->pivots == numpp) {
            dec_ctx->stage = FINALBACKWARD;
            finish_recovering_PP(dec_ctx);
        }
        return;
    }

    if (dec_ctx->stage == FORWARD && dec_ctx->pivots == dec_ctx->sc->snum) {
        // Attempt final forward substitution
        dec_ctx->stage = FINALFORWARD;
        dec_ctx->pivots = numpp - gensize;  // reset number of pivots before we verify the bottom rows
        // Transform the bottom wrap-around vectors to full length (numpp) and then process against their above band vectors
        // Copy the corresponding part of the message matrix as well
        GF_ELEMENT **ces = calloc(gensize, sizeof(GF_ELEMENT*));
        GF_ELEMENT **message = calloc(gensize, sizeof(GF_ELEMENT*));
        for (i=0; i<gensize; i++) {
            ces[i] = calloc(numpp, sizeof(GF_ELEMENT));
            if (ces[i] == NULL)
                fprintf(stderr, "%s: calloc ces failed\n", fname);
            for (j=0; j<dec_ctx->row[numpp-gensize+i]->len; j++) {
                int index = (numpp - gensize + i + j) % numpp;
                ces[i][index] = dec_ctx->row[numpp-gensize+i]->elem[j];
            }
            message[i] = calloc(pktsize, sizeof(GF_ELEMENT));
            memcpy(message[i], dec_ctx->message[numpp-gensize+i], pktsize*sizeof(GF_ELEMENT));
            // free the rows in dec_ctx->row, and message
            free(dec_ctx->row[numpp-gensize+i]->elem);
            free(dec_ctx->row[numpp-gensize+i]);
            dec_ctx->row[numpp-gensize+i] = NULL;
            memset(dec_ctx->message[numpp-gensize+i], 0, sizeof(GF_ELEMENT)*pktsize);
        }
        // Process one by one against above band vectors
        for (i=0; i<gensize; i++) {
            // Process the full-length row against above band vectors
            for (k=0; k<numpp; k++) {
                if (ces[i][k] != 0) {
                    if (dec_ctx->row[k] != NULL) {
                        assert(dec_ctx->row[k]->elem[0]);
                        quotient = galois_divide(ces[i][k], dec_ctx->row[k]->elem[0]);
                        galois_multiply_add_region(&(ces[i][k]), dec_ctx->row[k]->elem, quotient, dec_ctx->row[k]->len);
                        galois_multiply_add_region(message[i], dec_ctx->message[k], quotient, pktsize);
                        dec_ctx->operations += 1 + dec_ctx->row[k]->len + pktsize;
                    } else {
                        // a valid pivot found, store it back to decoding matrix
                        dec_ctx->row[k] = (struct row_vector*) malloc(sizeof(struct row_vector));
                        if (dec_ctx->row[k] == NULL)
                            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, k);
                        int len = numpp - k;
                        dec_ctx->row[k]->len = len;
                        dec_ctx->row[k]->elem = (GF_ELEMENT *) calloc(len, sizeof(GF_ELEMENT));
                        if (dec_ctx->row[k]->elem == NULL)
                            fprintf(stderr, "%s: calloc dec_ctx->row[%d]->elem failed\n", fname, k);
                        memcpy(dec_ctx->row[k]->elem, &(ces[i][k]), len*sizeof(GF_ELEMENT));
                        assert(dec_ctx->row[k]->elem[0]);
                        memcpy(dec_ctx->message[k], message[i], pktsize*sizeof(GF_ELEMENT));
                        dec_ctx->pivots += 1;
                        break;
                    }
                }
            }
            free(ces[i]);  // free the processed full-length vector
            free(message[i]);
        }
        // free ces, and message
        free(ces);
        free(message);  // this is just an array of pointers
    }
    if (dec_ctx->pivots == numpp) {
        dec_ctx->stage = FINALBACKWARD;
        finish_recovering_PP(dec_ctx);
    }
    return;
}

/**
 * Finish perpetual code decoding
 * This routine converts decoding matrix from upper triangular
 * form to diagonal.
 */
static void finish_recovering_PP(struct decoding_context_PP *dec_ctx)
{
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum;
    int i, j;
    int len;
    GF_ELEMENT quotient;
    // eliminate all nonzeros of the upper-triangular, examining columns from right to left
    for (i=numpp-1; i>=0; i--) {
        // examine rows from top to bottom
        // FIXME: j=0 is overkill, maybe j=max(0, i-gensize)?
        for (j=0; j<i; j++) {
            len = dec_ctx->row[j]->len;
            if (j+len <= i || dec_ctx->row[j]->elem[i-j] == 0)
                continue;
            assert(dec_ctx->row[i]->elem[0]);
            quotient = galois_divide(dec_ctx->row[j]->elem[i-j], dec_ctx->row[i]->elem[0]);
            galois_multiply_add_region(dec_ctx->message[j], dec_ctx->message[i], quotient, pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->row[j]->elem[i-j] = 0;
        }
        /* convert diagonal to 1*/
        if (dec_ctx->row[i]->elem[0] != 1) {
            galois_multiply_region(dec_ctx->message[i], galois_divide(1, dec_ctx->row[i]->elem[0]), pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->row[i]->elem[0] = 1;
        }
        /* save decoded packet */
        dec_ctx->sc->pp[i] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->sc->pp[i], dec_ctx->message[i], pktsize*sizeof(GF_ELEMENT));
    }
    dec_ctx->finished = 1;
}

void free_dec_context_PP(struct decoding_context_PP *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
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
    if (dec_ctx->sc != NULL)
        snc_free_enc_context(dec_ctx->sc);
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
long save_dec_context_PP(struct decoding_context_PP *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = PP_DECODER;
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
    filesize += fwrite(&dec_ctx->stage, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->pivots, sizeof(int), 1, fp);

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
    //performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

/*
 * Restore a decoding context from a file
 */
struct decoding_context_PP *restore_dec_context_PP(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    // Create a fresh decoding context
    struct decoding_context_PP *dec_ctx = create_dec_context_PP(&sp);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_CBD failed\n");
        return NULL;
    }
    // Restore decoding context from file
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field

    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->stage, sizeof(int), 1, fp);
    fread(&dec_ctx->pivots, sizeof(int), 1, fp);
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
                free_dec_context_PP(dec_ctx);
                return NULL;
            }
            dec_ctx->row[i]->len = rowlen;
            dec_ctx->row[i]->elem = (GF_ELEMENT *) calloc(rowlen, sizeof(GF_ELEMENT));
            if (dec_ctx->row[i]->elem == NULL) {
                free_dec_context_PP(dec_ctx);
                return NULL;
            }
            fread(dec_ctx->row[i]->elem, sizeof(GF_ELEMENT), rowlen, fp);
            fread(dec_ctx->message[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    //performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
