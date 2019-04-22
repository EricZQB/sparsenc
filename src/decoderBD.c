/*-----------------------decoderBD.c----------------------
 * Implementation of regular band decoder. It employs pivoting
 * to jointly decode band GNC code and its precode.
 *------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderBD.h"
static int partially_diag_decoding_matrix(struct decoding_context_BD *dec_ctx);
static int apply_parity_check_matrix(struct decoding_context_BD *dec_ctx);
static void finish_recovering_BD(struct decoding_context_BD *dec_ctx);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ctoo_r, int **ctoo_c, int *inactives);
extern long pivot_matrix_tworound(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ctoo_r, int **ctoo_c, int *inactives);

// create decoding context for band decoder
struct decoding_context_BD *create_dec_context_BD(struct snc_parameters *sp)
{
    static char fname[] = "snc_create_dec_context_BD";
    int i, j, k;

    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    if (sp->type != BAND_SNC) {
        fprintf(stderr, "Band decoder only applies to band GNC code.\n");
        return NULL;
    }

    struct decoding_context_BD *dec_ctx;
    if ((dec_ctx = malloc(sizeof(struct decoding_context_BD))) == NULL) {
        fprintf(stderr, "malloc decoding_context_BD failed\n");
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
    dec_ctx->inactivated  = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    dec_ctx->coefficient = calloc(numpp, sizeof(GF_ELEMENT*));
    if (dec_ctx->coefficient == NULL)
        goto AllocError;
    dec_ctx->message     = calloc(numpp, sizeof(GF_ELEMENT*));
    if (dec_ctx->message == NULL)
        goto AllocError;
    for (i=0; i<numpp; i++) {
        dec_ctx->coefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
        if (dec_ctx->coefficient[i] == NULL)
            goto AllocError;
        dec_ctx->message[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
        if (dec_ctx->message[i] == NULL)
            goto AllocError;
    }

    dec_ctx->ctoo_r = malloc(sizeof(int) * numpp);
    if (dec_ctx->ctoo_r == NULL)
        goto AllocError;
    dec_ctx->ctoo_c = malloc(sizeof(int) * numpp);
    if (dec_ctx->ctoo_c == NULL)
        goto AllocError;

    dec_ctx->overhead     = 0;
    dec_ctx->overheads = calloc(dec_ctx->sc->gnum, sizeof(int));
    if (dec_ctx->overheads == NULL)
        goto AllocError;
    dec_ctx->operations   = 0;
    return dec_ctx;

AllocError:
    free_dec_context_BD(dec_ctx);
    dec_ctx = NULL;
    return NULL;
}


void process_packet_BD(struct decoding_context_BD *dec_ctx, struct snc_packet *pkt)
{
    static char fname[] = "process_packet_BD";
    dec_ctx->overhead += 1;
    dec_ctx->overheads[pkt->gid] += 1;
    int i, j, k;
    int pivot;
    int pivotfound = 0;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // translate GNC encoding vector to full length
    GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));

    if (dec_ctx->de_precode == 0) {
        // Before precode's check matrix was applied.

        // Pay attention to systematic packets. It is easy to handle in BD decoder. Let's 
        // just transform it to full-length singleton vector.
        if (pkt->gid == -1 && pkt->ucid == -1) {
            fprintf(stderr, "%s: pkt's gid is -1 but ucid is not valid\n", fname);
            return;
        } else if (pkt->gid == -1 && pkt->ucid >= 0) {
            ces[pkt->ucid] = 1;
        } else {
            for (i=0; i<gensize; i++) {
                int index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
                if (dec_ctx->sc->params.gfpower==1) {
                    ces[index] = get_bit_in_array(pkt->coes, i);
                } else if (dec_ctx->sc->params.gfpower==8) {
                    ces[index] = pkt->coes[i];
                } else {
                    ces[index] = read_bits_from_byte_array(pkt->coes, dec_ctx->sc->params.size_g, dec_ctx->sc->params.gfpower, i);
                }
            }
        }
        for (i=0; i<numpp; i++) {
            if (ces[i] != 0) {
                if (dec_ctx->coefficient[i][i] != 0) {
                    quotient = galois_divide(ces[i], dec_ctx->coefficient[i][i]);
                    dec_ctx->operations += 1;
                    int band_width = numpp-i > gensize ? gensize : numpp-i;
                    galois_multiply_add_region(ces+i, &(dec_ctx->coefficient[i][i]), quotient, band_width);
                    dec_ctx->operations += band_width;
                    galois_multiply_add_region(pkt->syms, dec_ctx->message[i], quotient, pktsize);
                    dec_ctx->operations += pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
    } else {
        // Parity-check matrix has been applied, and therefore the decoding matrix has been pivoted and re-ordered
 
        // Pay attention to systematic packets. It is easy to handle in BD decoder. Let's 
        // just transform it to full-length singleton vector. Note that, however, it's very
        // unlikely that systematic packets would be received after parity-check matrix had
        // been applied.
        if (pkt->gid == -1 && pkt->ucid == -1) {
            fprintf(stderr, "%s: pkt's gid is -1 but ucid is not valid\n", fname);
            return;
        } else if (pkt->gid == -1 && pkt->ucid >= 0) {
            ces[pkt->ucid] = 1;
        } else {
            for (i=0; i<gensize; i++) {
                int orig_index = dec_ctx->sc->gene[pkt->gid]->pktid[i];
                if (dec_ctx->sc->params.gfpower == 1) {
                    ces[orig_index] = get_bit_in_array(pkt->coes, i);
                } else if (dec_ctx->sc->params.gfpower == 8){
                    ces[orig_index] = pkt->coes[i];
                } else {
                    ces[orig_index] = read_bits_from_byte_array(pkt->coes, dec_ctx->sc->params.size_g, dec_ctx->sc->params.gfpower, i);
                }
            }
        }
        for (i=0; i<numpp; i++) {
            if (ces[dec_ctx->ctoo_c[i]] != 0) {
                if (dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] != 0) {
                    quotient = galois_divide(ces[dec_ctx->ctoo_c[i]], dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]]);
                    dec_ctx->operations += 1;
                    for (j=i; j<numpp; j++) {
                        ces[dec_ctx->ctoo_c[j]] = galois_add(ces[dec_ctx->ctoo_c[j]], galois_multiply(dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]], quotient));
                    }
                    dec_ctx->operations += (numpp - i);
                    galois_multiply_add_region(pkt->syms, dec_ctx->message[dec_ctx->ctoo_r[i]], quotient, pktsize);
                    dec_ctx->operations += pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
    }

    //printf("Operations consumed: %lld\n", dec_ctx->operations);
     if (pivotfound == 1 && !dec_ctx->de_precode) {
        memcpy(dec_ctx->coefficient[pivot], ces, numpp*sizeof(GF_ELEMENT));
        memcpy(dec_ctx->message[pivot], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
        dec_ctx->DoF += 1;
     }

    if (pivotfound == 1 && dec_ctx->de_precode) {
        memcpy(dec_ctx->coefficient[dec_ctx->ctoo_r[pivot]], ces, numpp*sizeof(GF_ELEMENT));
        memcpy(dec_ctx->message[dec_ctx->ctoo_r[pivot]], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
        dec_ctx->DoF += 1;
    }
    // If the number of received DoF is equal to NUM_SRC, apply the parity-check matrix.
    // The messages corresponding to rows of parity-check matrix are all-zero.
    if (dec_ctx->de_precode == 0 && dec_ctx->DoF == dec_ctx->sc->snum) {
        if (get_loglevel() == TRACE)
            printf("Start to apply the parity-check matrix...\n");
        int allzeros = partially_diag_decoding_matrix(dec_ctx);
        if (get_loglevel() == TRACE)
            printf("%d all-zero rows when partially diagonalizing the decoding matrix.\n", allzeros);
        int missing_DoF = apply_parity_check_matrix(dec_ctx);
        if (get_loglevel() == TRACE)
            printf("After applying the parity-check matrix, %d DoF are missing.\n", missing_DoF);
        dec_ctx->DoF = numpp - missing_DoF;
        dec_ctx->de_precode = 1;
    }

    if (dec_ctx->DoF == dec_ctx->sc->snum + dec_ctx->sc->cnum) {
        finish_recovering_BD(dec_ctx);
    }

    free(ces);
    ces = NULL;
}

/*
 * Partially diagonalize the upper-trianguler decoding matrix,
 * i.e., remove nonzero elements above nonzero diagonal elements:
 *
 *   |x x x x       |      | x o x o       |
 *   |  x x x x     |      |   x x o o     |
 *   |    o o o o   |      |     o o o o   |
 *   |      x x x x |  --> |       x o x o |
 *   |        x x x |      |         x x o |
 *   |          o o |      |           o o |
 *   |            x |      |             x |
 *
 */
static int partially_diag_decoding_matrix(struct decoding_context_BD *dec_ctx)
{
    int         i, j, l;
    int         nonzero_rows = 0;
    GF_ELEMENT  quotient;
    long long   operations = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    int zero_size   = dec_ctx->sc->cnum + 5;
    int *zeropivots = (int *) malloc( sizeof(int) * zero_size );    // store indices of columns where the diagonal element is zero
    int zero_p      = 0;                                        // indicate how many zero pivots have been identified

    for (j=numpp-1; j>=0; j--) {
        if (dec_ctx->coefficient[j][j] == 0) {
            if (zero_p == zero_size) {
                zeropivots = (int *) realloc(zeropivots, sizeof(int) * (zero_size + 5));
                zero_size += 5;
            }
            zeropivots[zero_p++] = j;
            continue;
        } else {
            nonzero_rows += 1;
            int start_row = j-gensize > 0 ? j-gensize : 0;      // the upper triangular form is also in banded form, so no need to go through all rows
            for (i=start_row; i<j; i++) {
                if (dec_ctx->coefficient[i][j] == 0)
                    continue;

                quotient = galois_divide(dec_ctx->coefficient[i][j], dec_ctx->coefficient[j][j]);
                operations += 1;
                dec_ctx->coefficient[i][j] = 0;         // eliminiate the element
                // Important: corresponding operations on behind columns whose diagonal elements are zeros
                for (int z=0; z<zero_p; z++) {
                    l = zeropivots[z];
                    if (dec_ctx->coefficient[j][l] != 0) {
                        dec_ctx->coefficient[i][l] = galois_add(dec_ctx->coefficient[i][l], galois_multiply(dec_ctx->coefficient[j][l], quotient));
                        operations += 1;
                    }
                }
                // correspoding operations on the message matrix
                galois_multiply_add_region(dec_ctx->message[i], dec_ctx->message[j], quotient, pktsize);
                operations += pktsize;
            }
        }
    }
    free(zeropivots);
    dec_ctx->operations += operations;
    return (numpp-nonzero_rows);
}

// Apply the parity-check matrix to the decoding matrix; pivot, re-order and try to jointly decode
static int apply_parity_check_matrix(struct decoding_context_BD *dec_ctx)
{
    int i, j, k;
    int num_of_new_DoF = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
    int p = 0;                      // index pointer to the parity-check vector that is to be copyed
    for (i=0; i<numpp; i++) {
        if (dec_ctx->coefficient[i][i] == 0) {
            /* Set the coding vector according to parity-check bits */
            NBR_node *varnode = dec_ctx->sc->graph->l_nbrs_of_r[p]->first;
            while (varnode != NULL) {
                dec_ctx->coefficient[i][varnode->data] = varnode->ce;
                varnode = varnode->next;
            }
            dec_ctx->coefficient[i][dec_ctx->sc->snum+p] = 1;
            p++;
            memset(dec_ctx->message[i], 0, sizeof(GF_ELEMENT)*pktsize);         // parity-check vector corresponds to all-zero message
        }
    }

    // 2, Pivot and re-order matrices
    if (getenv("SNC_BD_ONEROUND") != NULL 
         && atoi(getenv("SNC_BD_ONEROUND")) == 1) {
        dec_ctx->operations += pivot_matrix_oneround(numpp, numpp, pktsize, dec_ctx->coefficient, dec_ctx->message, &dec_ctx->ctoo_r, &dec_ctx->ctoo_c, &(dec_ctx->inactivated));
    } else {
        dec_ctx->operations += pivot_matrix_tworound(numpp, numpp, pktsize, dec_ctx->coefficient, dec_ctx->message, &dec_ctx->ctoo_r, &dec_ctx->ctoo_c, &(dec_ctx->inactivated));
    }

    /* Count available innovative rows */
    int missing_DoF = 0;
    for (i=0; i<numpp; i++) {
        if (dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] == 0)
            missing_DoF++;
    }
    if (get_loglevel() == TRACE)
        printf("Missing %d DoF after applying parity-check matrix.\n", missing_DoF);
    return missing_DoF;
}


// recover decoded packets after NUM_SRC DoF has been received
static void finish_recovering_BD(struct decoding_context_BD *dec_ctx)
{
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    int i, j, k;
    GF_ELEMENT quotient;
    long long bs_ops = 0;
    // Backard substitution from right-most col to the left
    for (j=numpp-1; j>=0; j--) {
        for (i=0; i<j; i++) {
            if (dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]] != 0) {
                quotient = galois_divide(dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]], dec_ctx->coefficient[dec_ctx->ctoo_r[j]][dec_ctx->ctoo_c[j]]);
                galois_multiply_add_region(dec_ctx->message[dec_ctx->ctoo_r[i]], dec_ctx->message[dec_ctx->ctoo_r[j]], quotient, pktsize);
                dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]] = 0;
                bs_ops += 1 + pktsize;
            }
        }
    }
    // Convert all diagonal element to 1
    for (i=0; i<numpp; i++) {
        if (dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] != 1)
            galois_multiply_region(dec_ctx->message[dec_ctx->ctoo_r[i]], galois_divide(1, dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]]), pktsize);
        dec_ctx->coefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] = 1;
        int pktid = dec_ctx->ctoo_c[i];
        dec_ctx->sc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dec_ctx->sc->pp[pktid], dec_ctx->message[dec_ctx->ctoo_r[i]], pktsize*sizeof(GF_ELEMENT));
    }
    dec_ctx->operations += bs_ops;
    dec_ctx->finished = 1;
}

void free_dec_context_BD(struct decoding_context_BD *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    if (dec_ctx->coefficient != NULL) {
        for (int i=dec_ctx->sc->snum+dec_ctx->sc->cnum-1; i>=0; i--) {
            if (dec_ctx->coefficient[i] != NULL)
                free(dec_ctx->coefficient[i]);
        }
        free(dec_ctx->coefficient);
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
    if (dec_ctx->ctoo_r != NULL)
        free(dec_ctx->ctoo_r);
    if (dec_ctx->ctoo_c != NULL)
        free(dec_ctx->ctoo_c);
    if (dec_ctx->overheads != NULL)
        free(dec_ctx->overheads);
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
long save_dec_context_BD(struct decoding_context_BD *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = BD_DECODER;
    FILE *fp;
    if ((fp = fopen(filepath, "w")) == NULL) {
        fprintf(stderr, "Cannot open %s to save decoding context\n", filepath);
        return (-1);
    }
    int i, j;
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    // Write snc params
    filesize += fwrite(&dec_ctx->sc->params, sizeof(struct snc_parameters), 1, fp);
    // Write decoder type
    filesize += fwrite(&(d_type), sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->DoF, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->de_precode, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->inactivated, sizeof(int), 1, fp);
    // Save running matrices
    if (dec_ctx->de_precode == 0) {
        // Parity-check packets are not applied yet
        int count = 0;
        int len;
        i = 0;
        while (count != dec_ctx->DoF) {
            if (dec_ctx->coefficient[i][i] != 0) {
                filesize += fwrite(&i, sizeof(int), 1, fp);
                len = numpp -i < gensize ? numpp - i : gensize;
                filesize += fwrite(&len, sizeof(int), 1, fp);
                filesize += fwrite(&(dec_ctx->coefficient[i][i]), sizeof(GF_ELEMENT), len, fp);
                filesize += fwrite(dec_ctx->message[i], sizeof(GF_ELEMENT), pktsize, fp);
                count++;
            }
            i++;
        }
    } else {
        // Save the pivoted decoding matrix
        // Fixme: save space by making use of the sparse structure
        for (i=0; i<numpp; i++) {
            filesize += fwrite(dec_ctx->coefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            filesize += fwrite(dec_ctx->message[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
        filesize += fwrite(dec_ctx->ctoo_r, sizeof(int), numpp, fp);
        filesize += fwrite(dec_ctx->ctoo_c, sizeof(int), numpp, fp);
    }
    // Save performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(dec_ctx->overheads, sizeof(int), dec_ctx->sc->gnum, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}


struct decoding_context_BD *restore_dec_context_BD(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field
    // Create a fresh decoding context
    struct decoding_context_BD *dec_ctx = create_dec_context_BD(&sp);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_GG failed\n");
        return NULL;
    }
    // Restore decoding context from file
    int i, j;
    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->DoF, sizeof(int), 1, fp);
    fread(&dec_ctx->de_precode, sizeof(int), 1, fp);
    fread(&dec_ctx->inactivated, sizeof(int), 1, fp);
    // Restore running matrices
    if (dec_ctx->de_precode == 0) {
        int count = 0;
        int pivot, len;
        while (count != dec_ctx->DoF) {
           fread(&pivot, sizeof(int), 1, fp);
           fread(&len, sizeof(int), 1, fp);
           fread(&(dec_ctx->coefficient[pivot][pivot]), sizeof(GF_ELEMENT), len, fp);
           fread(dec_ctx->message[pivot], sizeof(GF_ELEMENT), sp.size_p, fp);
           count++;
        }
    } else {
        int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;
        for (i=0; i<numpp; i++) {
            fread(dec_ctx->coefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            fread(dec_ctx->message[i], sizeof(GF_ELEMENT), numpp, fp);
        }
        fread(dec_ctx->ctoo_r, sizeof(int), numpp, fp);
        fread(dec_ctx->ctoo_c, sizeof(int), numpp, fp);
    }
    // Restore performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(dec_ctx->overheads, sizeof(int), dec_ctx->sc->gnum, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
