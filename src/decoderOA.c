/*--------------------------sncOADecoder.c----------------------
 * Implementation of overlap-aware (OA) decoder.
 *-------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderOA.h"

struct running_matrix
{
    struct row_vector **row;
    GF_ELEMENT **message;
};
/*
 * Transform coefficient matrix of each generation to row echelon form (REF).
 * The REF has the property that if a diagonal element is nonzero, all elements
 * above it are reduced to zero.
 */
static long running_matrix_to_REF(struct decoding_context_OA *dec_ctx);

/*
 * Construct global decoding matrix (GDM) from REFs of each generation. Two times of
 * pivoting are done to transform GDM to a specific form where the bottom-right corner is
 * dense and the above half is partially diagonal:
 * -                 -
 * | x 0 0 0 0 x x x |
 * | 0 x 0 0 0 x x x |
 * | 0 0 x 0 0 x x x |
 * | 0 0 0 x 0 x x x |
 * | 0 0 0 0 x x x x |
 * | 0 0 0 0 0 x x x |
 * | 0 0 0 0 0 x x x |
 * | 0 0 0 0 0 x x x |
 * -                -
 */
static void construct_GDM(struct decoding_context_OA *dec_ctx);

/*
 * Fully transform GDM to identity matrix to finish decoding.
 */
static void diagonalize_GDM(struct decoding_context_OA *dec_ctx);

/* Free running matrix */
static void free_running_matrix(struct running_matrix *mat, int rows);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ctoo_r, int **ctoo_c, int *inactives);
extern long pivot_matrix_tworound(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ctoo_r, int **ctoo_c, int *inactives);

static clock_t l_proc_start = 0;   // local processing start

/*
 * snc_create_dec_context_OA
 * Create context for overlap-aware (OA) decoding
 *  aoh - allowed overhead >=0
 */
struct decoding_context_OA *create_dec_context_OA(struct snc_parameters *sp, int aoh)
{
    static char fname[] = "snc_create_dec_context_OA";
    int i, j, k;

    struct decoding_context_OA *dec_ctx;
    if ((dec_ctx = malloc(sizeof(struct decoding_context_OA))) == NULL)
        return NULL;

    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }

    dec_ctx->sc = sc;

    dec_ctx->aoh        = aoh;
    dec_ctx->finished   = 0;
    dec_ctx->OA_ready   = 0;
    dec_ctx->local_DoF  = 0;
    dec_ctx->global_DoF = 0;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // Allocate matrices for per-generation decoding
    dec_ctx->Matrices = calloc(dec_ctx->sc->gnum, sizeof(struct running_matrix*));
    if (dec_ctx->Matrices == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->Matrices\n", fname);
        goto AllocError;
    }
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        dec_ctx->Matrices[i] = calloc(1, sizeof(struct running_matrix));
        if (dec_ctx->Matrices[i] == NULL) {
            fprintf(stderr, "%s: malloc dec_ctx->Matrices[%d]\n", fname, i);
            goto AllocError;
        }
        // Allocate coefficient and message matrices in running_matrix
        // coefficeint: size_g x size_g
        // message:     size_g x size_p
        // Dim-1) Pointers to each row
        dec_ctx->Matrices[i]->row = calloc(gensize, sizeof(struct row_vector *));
        if (dec_ctx->Matrices[i]->row == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->row\n", fname, i);
            goto AllocError;
        }
        dec_ctx->Matrices[i]->message = calloc(gensize, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->message == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage\n", fname, i);
            goto AllocError;
        }
    }

    /*
     * We don't allocate memory for global decoding (ie GDM) here. We only allocate
     * when OA ready. This avoids occupying a big amount of memory for a long time.
     */

    // performance indices
    dec_ctx->operations = 0;
    dec_ctx->overhead   = 0;
    dec_ctx->ops1 = dec_ctx->ops2 = dec_ctx->ops3 = dec_ctx->ops4 = 0;
    return dec_ctx;

AllocError:
    free_dec_context_OA(dec_ctx);
    dec_ctx = NULL;
    return NULL;
}

void process_packet_OA(struct decoding_context_OA *dec_ctx, struct snc_packet *pkt)
{
    static char fname[] = "process_packet_OA";
    if (l_proc_start == 0) {
        l_proc_start = clock();  // initialize timer
    }
    dec_ctx->overhead += 1;

    int i, j, k;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // start processing
    int gid = pkt->gid;
    int pivotfound = 0;
    int pivot;

    /*
     * If decoder is not OA ready, process the packet within the generation.
     */
    if (dec_ctx->OA_ready != 1) {
        GF_ELEMENT *pkt_coes = calloc(gensize, sizeof(GF_ELEMENT));

        // Need to pay extra attention to systematic packets
        if (pkt->gid == -1 && pkt->ucid == -1) {
            fprintf(stderr, "%s: pkt's gid is -1 but ucid is not valid\n", fname);
        } else if (pkt->gid == -1 && pkt->ucid >= 0) {
            // systematic packet
            gid = pkt->ucid / dec_ctx->sc->params.size_b;   // process it as a packet belonging to its base subgeneration;
            int local_pos = pkt->ucid % dec_ctx->sc->params.size_b;
            pkt_coes[local_pos] = 1;
        } else {
            // normal GNC packets
            if (dec_ctx->sc->params.bnc) {
                for (i=0; i<gensize; i++)
                    pkt_coes[i] = get_bit_in_array(pkt->coes, i);
            } else {
                memcpy(pkt_coes, pkt->coes, gensize*sizeof(GF_ELEMENT));
            }
        }
        // Translate the encoding vector to the sorted form as in the generation
        struct running_matrix *matrix = dec_ctx->Matrices[gid];
        for (i=0; i<gensize; i++) {
            if (pkt_coes[i] != 0) {
                if (matrix->row[i] != NULL) {
                    quotient = galois_divide(pkt_coes[i], matrix->row[i]->elem[0]);
                    galois_multiply_add_region(&(pkt_coes[i]), matrix->row[i]->elem, quotient, matrix->row[i]->len);
                    galois_multiply_add_region(pkt->syms, matrix->message[i], quotient, pktsize);
                    dec_ctx->operations += 1 + matrix->row[i]->len + pktsize;
                    dec_ctx->ops1 += 1 + matrix->row[i]->len + pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
        // cache as normal GNC packet
        if (pivotfound == 1) {
            matrix->row[pivot] = malloc(sizeof(struct row_vector));
            matrix->row[pivot]->len = gensize - pivot;
            matrix->row[pivot]->elem = malloc(sizeof(GF_ELEMENT) * matrix->row[pivot]->len);
            memcpy(matrix->row[pivot]->elem, &(pkt_coes[pivot]), sizeof(GF_ELEMENT)*matrix->row[pivot]->len);
            matrix->message[pivot] = malloc(sizeof(GF_ELEMENT) * pktsize);
            memcpy(matrix->message[pivot], pkt->syms, pktsize*sizeof(GF_ELEMENT));
            dec_ctx->local_DoF += 1;
        }
        free(pkt_coes);

        if (dec_ctx->local_DoF >= (dec_ctx->sc->snum+dec_ctx->aoh)) {
            dec_ctx->OA_ready = 1;
            // When OA ready, convert LDMs to upper triangular form
            long ops = running_matrix_to_REF(dec_ctx);
            dec_ctx->operations += ops;
            dec_ctx->ops1 += ops;
            // Record time used between processing the first received packet and OA ready
            if (get_loglevel() == TRACE) {
                printf("Local processing took %.6f seconds\n", ((double) (clock()-l_proc_start))/CLOCKS_PER_SEC);
            }
            // Combine LDMs to GDM and apply inactivation pivoting
            clock_t start, stop;
            start = clock();
            construct_GDM(dec_ctx);
            stop = clock();
            if (get_loglevel() == TRACE) {
                printf("Construct GDM took %.6f seconds\n", ((double) (stop-start))/CLOCKS_PER_SEC);
            }

            // If numpp innovative packets are received, recover all
            // source packets from JMBcoeffcient and JMBmessage
            if (dec_ctx->global_DoF == numpp) {
                start = clock();
                diagonalize_GDM(dec_ctx);
                stop = clock();
                if (get_loglevel() == TRACE) {
                    printf("Diagonalize GDM took %.6f seconds\n", ((double) (stop-start))/CLOCKS_PER_SEC);
                }
            }
        }
    } else {
        /*
         * If decoder is OA ready, process the packet against the global matrix directly.
         * Since the decoder is OA ready, translate the local encoded vector (LEV)
         * to global encoding vector (GEV). Since the GDM was probably pivoted, need
         * to transform the GEV according to the pivoting order.
         */
        GF_ELEMENT *re_ordered = calloc(numpp, sizeof(GF_ELEMENT));
        for (i=0; i<gensize; i++) {
            /* obtain index position of pktid in the full-length vector */
            int curr_pos = dec_ctx->sc->gene[gid]->pktid[i];
            if (dec_ctx->sc->params.bnc) {
                re_ordered[curr_pos] = get_bit_in_array(pkt->coes, i);
            } else {
                re_ordered[curr_pos] = pkt->coes[i];
            }
        }

        /*
         * Process the reordered GEV against GDM
         */
        pivot = -1;
        for (int m=0; m<numpp; m++) {
            if (re_ordered[dec_ctx->ctoo_c[m]] != 0) {
                if (dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[m]][dec_ctx->ctoo_c[m]] != 0) {
                    // mask the encoding vector and message over the JMB decoding matrix
                    GF_ELEMENT quotient = galois_divide(re_ordered[dec_ctx->ctoo_c[m]], dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[m]][dec_ctx->ctoo_c[m]]);
                    dec_ctx->operations += 1;
                    dec_ctx->ops3 += 1;
                    for (j=m; j<numpp; j++) {
                        re_ordered[dec_ctx->ctoo_c[j]] = galois_add(re_ordered[dec_ctx->ctoo_c[j]], galois_multiply(dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[m]][dec_ctx->ctoo_c[j]], quotient));
                    }
                    galois_multiply_add_region(pkt->syms, dec_ctx->JMBmessage[dec_ctx->ctoo_r[m]], quotient, pktsize);
                    dec_ctx->operations += pktsize;
                    dec_ctx->ops3 += pktsize;
                } else {
                    pivotfound = 1;
                    pivot = m;
                    break;
                }
            }
        }

        if (pivotfound == 1) {
            memcpy(dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[pivot]], re_ordered, numpp*sizeof(GF_ELEMENT));
            memcpy(dec_ctx->JMBmessage[dec_ctx->ctoo_r[pivot]], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
            dec_ctx->global_DoF += 1;

            if (dec_ctx->global_DoF == numpp) {
                // recover all source from JMBcoeffcient & JMBmessage matrix
                clock_t start, stop;
                start = clock();
                diagonalize_GDM(dec_ctx);
                stop = clock();
                if (get_loglevel() == TRACE) {
                    printf("Diagonalize GDM took %.6f seconds\n", ((double) (stop-start))/CLOCKS_PER_SEC);
                }
            }
        }
        free(re_ordered);
    }
    if (dec_ctx->finished && get_loglevel() == TRACE) {
        printf("OA splitted operations: %.2f %.2f %.2f %.2f\n",
                (double) dec_ctx->ops1/dec_ctx->sc->snum/pktsize,
                (double) dec_ctx->ops2/dec_ctx->sc->snum/pktsize,
                (double) dec_ctx->ops3/dec_ctx->sc->snum/pktsize,
                (double) dec_ctx->ops4/dec_ctx->sc->snum/pktsize);
    }
}

void free_dec_context_OA(struct decoding_context_OA *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    if (dec_ctx->sc != NULL)
        snc_free_enc_context(dec_ctx->sc);
    int i, j, k;
    if (dec_ctx->Matrices != NULL) {
        for (i=0; i<dec_ctx->sc->gnum; i++){
            // Free each decoding matrix
            if (dec_ctx->Matrices[i] != NULL)
                free_running_matrix(dec_ctx->Matrices[i], dec_ctx->sc->params.size_g);
            dec_ctx->Matrices[i] = NULL;
        }
        free(dec_ctx->Matrices);
        dec_ctx->Matrices = NULL;
    }
    if (dec_ctx->JMBcoefficient != NULL) {
        for (j=0; j<dec_ctx->sc->snum+dec_ctx->sc->cnum+dec_ctx->aoh; j++) {
            if (dec_ctx->JMBcoefficient[j] != NULL)
                free(dec_ctx->JMBcoefficient[j]);
        }
        free(dec_ctx->JMBcoefficient);
    }
    if (dec_ctx->JMBmessage != NULL) {
        for (j=0; j<dec_ctx->sc->snum+dec_ctx->sc->cnum+dec_ctx->aoh; j++) {
            if (dec_ctx->JMBmessage[j] != NULL)
                free(dec_ctx->JMBmessage[j]);
        }
        free(dec_ctx->JMBmessage);
    }
    if (dec_ctx->ctoo_r != NULL)
        free(dec_ctx->ctoo_r);
    if (dec_ctx->ctoo_c != NULL)
        free(dec_ctx->ctoo_c);
    free(dec_ctx);
    dec_ctx = NULL;
    return;
}

static void free_running_matrix(struct running_matrix *mat, int rows)
{
    if (mat != NULL) {
        for (int i=0; i<rows; i++) {
            if (mat->row[i] != NULL) {
                if (mat->row[i]->elem != NULL)
                    free(mat->row[i]->elem);
                free(mat->row[i]);
                mat->row[i] = NULL;
            }
        }
        if (mat->message != NULL) {
            for (int i=0; i<rows; i++) {
                if (mat->message[i] != NULL)
                    free(mat->message[i]);
                    mat->message[i] = NULL;
            }
            free(mat->message);
            mat->message = NULL;
        }
        free(mat);
        return;
    }
}

static void diagonalize_GDM(struct decoding_context_OA *dec_ctx)
{
    static char fname[] = "finish_recovering_inactivation";
    int i, j, k;
    int pos;
    int pktid;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    // Recover inactivated packets
    if (get_loglevel() == TRACE) {
        printf("Finishing decoding...\n");
        printf("Recovering \"inactive\" packets...\n");
    }
    int ias = dec_ctx->inactives;
    GF_ELEMENT **ces_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
    GF_ELEMENT **msg_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++) {
        ces_submatrix[i] = calloc(ias, sizeof(GF_ELEMENT));
        msg_submatrix[i] = calloc(pktsize, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            ces_submatrix[i][j] = dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[numpp-ias+i]][dec_ctx->ctoo_c[numpp-ias+j]];
        memcpy(msg_submatrix[i], dec_ctx->JMBmessage[dec_ctx->ctoo_r[numpp-ias+i]], pktsize*sizeof(GF_ELEMENT));
    }

    /* Perform back substitution to reduce the "ias x ias" matrix to identity matrix */
    long long ops = back_substitute(ias, ias, pktsize, ces_submatrix, msg_submatrix);
    dec_ctx->operations += ops;
    dec_ctx->ops4 += ops;

    // Recover decoded overlapping packets
    for (i=0; i<ias; i++) {
        // get original pktid at column (numpp-ias+i0
        pktid = dec_ctx->ctoo_c[numpp-ias+i];
        // Construct decoded packets
        if ( (dec_ctx->sc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
            fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, pktid);
        memcpy(dec_ctx->sc->pp[pktid], msg_submatrix[i], sizeof(GF_ELEMENT)*pktsize);
    }
    // Copy back msg_submatrix
    for (i=0; i<ias; i++)
        memcpy(dec_ctx->JMBmessage[dec_ctx->ctoo_r[numpp-ias+i]], msg_submatrix[i], pktsize*sizeof(GF_ELEMENT));
    // free ces_submatrix, msg_submatrix
    for (i=0; i<ias; i++) {
        free(ces_submatrix[i]);
        free(msg_submatrix[i]);
    }
    free(ces_submatrix);
    free(msg_submatrix);


    // Recover active packets
    if (get_loglevel() == TRACE)
        printf("Recovering \"active\" packets...\n");
    GF_ELEMENT quotient;
    for (i=0; i<numpp-ias; i++) {
        /*
         * Clean up the inactive part of the upper half of GDM by
         * masking non-zero element aginst already decoded inactive packets
         *
         */
        for (j=numpp-ias; j<numpp; j++) {
            if (dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]] != 0) {
                quotient = dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]];
                pktid = dec_ctx->ctoo_c[j];
                galois_multiply_add_region(dec_ctx->JMBmessage[dec_ctx->ctoo_r[i]], dec_ctx->sc->pp[pktid], quotient, pktsize);
                dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[j]] = 0;
                dec_ctx->operations += pktsize;
                dec_ctx->ops4 += pktsize;
            }
        }

        // Convert diagonal elements of top-left part of T to 1
        quotient = dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]];
        if (quotient != 1) {
            galois_multiply_region(dec_ctx->JMBmessage[dec_ctx->ctoo_r[i]], galois_divide(1, quotient), pktsize);
            dec_ctx->operations += pktsize;
            dec_ctx->ops4 += pktsize;
            dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] = 1;
        }

        // Save the decoded packet
        pktid = dec_ctx->ctoo_c[i];
        if ( dec_ctx->sc->pp[pktid] != NULL )
            fprintf(stderr, "%s：warning: packet %d is already recovered.\n", fname, pktid);
        if ( (dec_ctx->sc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
            fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, pktid);
        memcpy(dec_ctx->sc->pp[pktid], dec_ctx->JMBmessage[dec_ctx->ctoo_r[i]], sizeof(GF_ELEMENT)*pktsize);
    }

    dec_ctx->finished = 1;
}

// Partially diagonalize all running matrices when the decoder is OA ready
// diagonalizes them as much as possible
static long running_matrix_to_REF(struct decoding_context_OA *dec_ctx)
{
    long long operations = 0;
    int i, j, k, l;
    GF_ELEMENT quotient;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;


    #pragma omp parallel for private(i)
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        struct running_matrix *matrix = dec_ctx->Matrices[i];

        // Partially diagonalize the LDM
        int consecutive = 1;
        for (k=gensize-1; k>=0; k--) {
            if (matrix->row[k] == NULL) {
                consecutive = 0;
                continue;
            }

            // eliminate elements above the nonzero diagonal elements
            for (l=0; l<k; l++) {
                if (matrix->row[l] == NULL || matrix->row[l]->elem[k-l] == 0)
                    continue;

                quotient = galois_divide(matrix->row[l]->elem[k-l], matrix->row[k]->elem[0]);
                operations += 1;
                matrix->row[l]->elem[k-l] = 0;
                // Note that columns behind the current column could be nonzero because of their zero diagonal.
                for (int m=k+1; m<gensize; m++) {
                    if (matrix->row[m] == NULL) {
                        matrix->row[l]->elem[m-l] = galois_add(matrix->row[l]->elem[m-l], galois_multiply(matrix->row[k]->elem[m-k], quotient));
                        operations += 1;
                    }
                }
                galois_multiply_add_region(matrix->message[l], matrix->message[k], quotient, pktsize);
                operations += pktsize;
            }
        }

        if (consecutive == 1 && get_loglevel() == TRACE)
            printf("Class %d is self-decodable.\n", i);
    }
    return operations;
}

static void construct_GDM(struct decoding_context_OA *dec_ctx)
{
    static char fname[] = "construct_GDM";
    int i, j, k;
    struct running_matrix *matrix;

    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;

    //Allocate GDM to snc_dec_context, apply precoding matrix
    dec_ctx->JMBcoefficient = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
    if (dec_ctx->JMBcoefficient == NULL)
        fprintf(stderr, "%s: calloc dec_ctx->JMBcoefficient\n", fname);
    dec_ctx->JMBmessage     = calloc(numpp+dec_ctx->aoh, sizeof(GF_ELEMENT*));
    if (dec_ctx->JMBmessage == NULL)
        fprintf(stderr, "%s: calloc dec_ctx->JMBmessage\n", fname);

    for (i=0; i<numpp+dec_ctx->aoh; i++) {
        dec_ctx->JMBcoefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
        dec_ctx->JMBmessage[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
    }
    dec_ctx->inactives   = 0;
    dec_ctx->ctoo_r = malloc(sizeof(int) * numpp);
    dec_ctx->ctoo_c = malloc(sizeof(int) * numpp);
    // Apply precoding matrix
    for (i=0; i<dec_ctx->sc->cnum; i++) {
        dec_ctx->JMBcoefficient[dec_ctx->sc->snum+dec_ctx->aoh+i][dec_ctx->sc->snum+i] = 1;

        NBR_node *variable_node = dec_ctx->sc->graph->l_nbrs_of_r[i]->first;        //ldpc_graph->nbrs_of_right[i];
        while (variable_node != NULL) {
            // 标记与该check packet连结的所有source packet node
            int src_pktid = variable_node->data;                        //variable_node->nb_index;
            dec_ctx->JMBcoefficient[dec_ctx->sc->snum+dec_ctx->aoh+i][src_pktid] = variable_node->ce;
            variable_node = variable_node->next;
        }
    }

    // Step 1, translate LEVs to GEV and move them to GDM
    GF_ELEMENT *global_ces = calloc(numpp, sizeof(GF_ELEMENT));
    int p_copy = 0;                             // 拷贝到JMBcofficient的行指针
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        matrix = dec_ctx->Matrices[i];
        for (j=0; j<gensize; j++) {
            if (matrix->row[j] == NULL)
                continue;                       // there is no local DoF here
            else {
                memset(global_ces, 0, numpp*sizeof(GF_ELEMENT));    /* Reset before reuse */
                for (k=j; k<gensize; k++)
                    global_ces[dec_ctx->sc->gene[i]->pktid[k]] = matrix->row[j]->elem[k-j];
                memcpy(dec_ctx->JMBcoefficient[p_copy], global_ces, numpp*sizeof(GF_ELEMENT));
                memcpy(dec_ctx->JMBmessage[p_copy], matrix->message[j], pktsize*sizeof(GF_ELEMENT));
                p_copy += 1;
            }
        }
    }
    free(global_ces);
    if (get_loglevel() == TRACE)
        printf("%d local DoFs are available, copied %d to GDM.\n", dec_ctx->local_DoF, p_copy);
    // Free up local matrices
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        free_running_matrix(dec_ctx->Matrices[i], dec_ctx->sc->params.size_g);
        dec_ctx->Matrices[i] == NULL;
    }
    free(dec_ctx->Matrices);
    dec_ctx->Matrices = NULL;

    /* Transform GDM to upper trianguler via pivoting */
    clock_t start_pivoting, stop_pivoting;
    start_pivoting = clock();
    long long ops;
    if (getenv("SNC_OA_ONEROUND") != NULL
         && atoi(getenv("SNC_OA_ONEROUND")) == 1) {
        ops = pivot_matrix_oneround(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, &dec_ctx->ctoo_r, &dec_ctx->ctoo_c, &(dec_ctx->inactives));
    } else {
        ops = pivot_matrix_tworound(numpp+dec_ctx->aoh, numpp, pktsize, dec_ctx->JMBcoefficient, dec_ctx->JMBmessage, &dec_ctx->ctoo_r, &dec_ctx->ctoo_c, &(dec_ctx->inactives));
    }
    stop_pivoting = clock();
    dec_ctx->operations += ops;
    dec_ctx->ops2 += ops;
    // Count available degree of freedom
    for (i=0; i<numpp; i++) {
        if (dec_ctx->JMBcoefficient[dec_ctx->ctoo_r[i]][dec_ctx->ctoo_c[i]] != 0)
            dec_ctx->global_DoF++;
    }
    double pivot_time = ((double) (stop_pivoting - start_pivoting)) / CLOCKS_PER_SEC;
    if (get_loglevel() == TRACE) {
        printf("Pivoting and forward substitution in total consumed %.6f seconds.\n", pivot_time);
        printf("A total of %d DoF have been received.\n", dec_ctx->global_DoF);
    }
}

/**
 * Save a decoding context to a file
 * Return values:
 *   On success: bytes written
 *   On error: -1
 */
long save_dec_context_OA(struct decoding_context_OA *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = OA_DECODER;
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
    filesize += fwrite(&dec_ctx->aoh, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->OA_ready, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->local_DoF, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->global_DoF, sizeof(int), 1, fp);
    // Save running matrices
    // FIXME: if OA_ready=1, Matrices have been freed.
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        for (j=0; j<gensize; j++) {
            int rowlen = dec_ctx->Matrices[i]->row[j] == NULL ? 0 : dec_ctx->Matrices[i]->row[j]->len;
            filesize += fwrite(&rowlen, sizeof(int), 1, fp);
            if (rowlen != 0) {
                filesize += fwrite(dec_ctx->Matrices[i]->row[j]->elem, sizeof(GF_ELEMENT), rowlen, fp);
                filesize += fwrite(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), pktsize, fp);
            }
        }
    }

    // Save GDM and its related bookkeeping information if OA ready.
    // Fixme: GDM is sparse, so we can save space by using compact representation
    // of the matrix. For example, save [row, col, value] for nonzero elements.
    if (dec_ctx->OA_ready == 1) {
        for (i=0; i<numpp + dec_ctx->aoh; i++) {
            filesize += fwrite(dec_ctx->JMBcoefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            filesize += fwrite(dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
        filesize += fwrite(dec_ctx->ctoo_r, sizeof(int), numpp, fp);
        filesize += fwrite(dec_ctx->ctoo_c, sizeof(int), numpp, fp);
        filesize += fwrite(&dec_ctx->inactives, sizeof(int), 1, fp);
    }
    // Save performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops1, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops2, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops3, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops4, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

struct decoding_context_OA *restore_dec_context_OA(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field
    int aoh;
    fread(&aoh, sizeof(int), 1, fp);
    // Create a fresh decoding context
    struct decoding_context_OA *dec_ctx = create_dec_context_OA(&sp, aoh);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_GG failed\n");
        return NULL;
    }
    // Restore decoding context from file
    int i, j;
    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->OA_ready, sizeof(int), 1, fp);
    fread(&dec_ctx->local_DoF, sizeof(int), 1, fp);
    fread(&dec_ctx->global_DoF, sizeof(int), 1, fp);
    // Restore running matrices
    // Note that running matrices' memory were already allocated in creating_dec_context
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        for (j=0; j<sp.size_g; j++) {
            int rowlen = 0;
            fread(&rowlen, sizeof(int), 1, fp);
            if (rowlen != 0) {
                dec_ctx->Matrices[i]->row[j] = malloc(sizeof(struct row_vector));
                dec_ctx->Matrices[i]->row[j]->len = rowlen;
                dec_ctx->Matrices[i]->row[j]->elem = (GF_ELEMENT *) malloc(rowlen * sizeof(GF_ELEMENT));
                fread(dec_ctx->Matrices[i]->row[j]->elem, sizeof(GF_ELEMENT), rowlen, fp);
                dec_ctx->Matrices[i]->message[j] = (GF_ELEMENT *) malloc(sp.size_p * sizeof(GF_ELEMENT));
                fread(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), sp.size_p, fp);
            }
        }
    }
    // Restore GDM and its related information if OA_ready
    int numpp = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    if (dec_ctx->OA_ready == 1) {
        dec_ctx->JMBcoefficient = calloc(numpp+aoh, sizeof(GF_ELEMENT*));
        dec_ctx->JMBmessage = calloc(numpp+aoh, sizeof(GF_ELEMENT*));
        for (i=0; i<numpp+aoh; i++) {
            dec_ctx->JMBcoefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
            fread(dec_ctx->JMBcoefficient[i], sizeof(GF_ELEMENT), numpp, fp);
            dec_ctx->JMBmessage[i] = calloc(sp.size_p, sizeof(GF_ELEMENT));
            fread(dec_ctx->JMBmessage[i], sizeof(GF_ELEMENT), sp.size_p, fp);
        }
        dec_ctx->ctoo_r = calloc(numpp, sizeof(int));
        fread(dec_ctx->ctoo_r, sizeof(int), numpp, fp);
        dec_ctx->ctoo_c = calloc(numpp, sizeof(int));
        fread(dec_ctx->ctoo_c, sizeof(int), numpp, fp);
        fread(&dec_ctx->inactives, sizeof(int), 1, fp);
    }
    // Restore performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops1, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops2, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops3, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops4, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}
