/*------------------------- decoderGG.c -----------------------
 * Implementation of generation-by-generation decoding.
 *-------------------------------------------------------------*/
#include "common.h"
#include "galois.h"
#include "decoderGG.h"

struct running_matrix {
    int DoF_miss;
    unsigned char *erased;   // bits indicating recovered packets
    struct row_vector **row;
    GF_ELEMENT **message;
};

static void decode_generation(struct decoding_context_GG *dec_ctx, int gid);
static void perform_iterative_decoding(struct decoding_context_GG *dec_ctx);
static void new_decoded_source_packet(struct decoding_context_GG *dec_ctx, int pkt_id);
static void new_decoded_check_packet(struct decoding_context_GG *dec_ctx, int pkt_id);
static void update_generations(struct decoding_context_GG *dec_ctx);
static long update_running_matrix(struct decoding_context_GG *dec_ctx, int gid, int sid, int index);
static int check_for_new_recoverables(struct decoding_context_GG *dec_ctx);
static int check_for_new_decodables(struct decoding_context_GG *dec_ctx);
static void free_running_matrix(struct running_matrix *matrix, int rows);

static ID_list **gene_nbr = NULL;    // lists of subgeneration neighbors of each packet

// Performance analysis use ONLY
static int g_decoded = 0; // decoded from subgen
static int c_decoded = 0; // recovered from precode decoding

// setup decoding context:
struct decoding_context_GG *create_dec_context_GG(struct snc_parameters *sp)
{
    static char fname[] = "snc_create_dec_context_GG";
    int i, j;

    struct decoding_context_GG *dec_ctx;
    if ((dec_ctx = malloc(sizeof(struct decoding_context_GG))) == NULL) {
        fprintf(stderr, "%s: malloc decoding context GG failed\n", fname);
        return NULL;
    }
    // GNC code context
    // Since this is decoding, we construct GNC context without data
    // sc->pp will be filled by decoded packets
    struct snc_context *sc;
    if ((sc = snc_create_enc_context(NULL, sp)) == NULL) {
        fprintf(stderr, "%s: create decoding context failed", fname);
        goto AllocError;
    }
    dec_ctx->sc = sc;

    // Since GG decoder frequently needs to find out which generations a packet belongs to, we
    // build the lists of subgeneration neighbors of each packet according to sc->gene. This list
    // avoids calling has_item() in all subgenerations.
    gene_nbr = build_subgen_nbr_list(sc);

    // memory areas needed for decoding
    dec_ctx->evolving_checks = calloc(dec_ctx->sc->cnum, sizeof(GF_ELEMENT *));
    if (dec_ctx->evolving_checks == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->evolving_checks\n", fname);
        goto AllocError;
    }
    dec_ctx->check_degrees = calloc(dec_ctx->sc->cnum, sizeof(int));
    if (dec_ctx->check_degrees == NULL) {
        fprintf(stderr, "%s: calloc dec_ctx->check_degrees\n", fname);
        goto AllocError;
    }
    for (i=0; i<dec_ctx->sc->cnum; i++) {
        NBR_node *nb = dec_ctx->sc->graph->l_nbrs_of_r[i]->first;
        while (nb != NULL) {
            dec_ctx->check_degrees[i]++;    // initial check degree of each check packet
            nb = nb->next;
        }
    }

    dec_ctx->finished  = 0;
    dec_ctx->decoded   = 0;
    dec_ctx->originals = 0;
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
        dec_ctx->Matrices[i]->DoF_miss = dec_ctx->sc->params.size_g;

        // Bits indicating column status
        int nflags = ALIGN(dec_ctx->sc->params.size_g, 8);  // Each byte has 8 bits
        dec_ctx->Matrices[i]->erased = calloc(nflags, sizeof(unsigned char));
        if (dec_ctx->Matrices[i]->erased == NULL) {
            fprintf(stderr, "%s: malloc dec_ctx->Matrices[%d]->erased\n", fname, i);
            goto AllocError;
        }

        // Allocate coefficient rows and message matrices in running_matrix
        // row: size_g rows, which are initialized to NULL
        // message:     size_g x size_p
        dec_ctx->Matrices[i]->row = (struct row_vector **) calloc(dec_ctx->sc->params.size_g, sizeof(struct row_vector *));
        if (dec_ctx->Matrices[i]->row == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->row failed\n", fname, i);
            goto AllocError;
        }
        dec_ctx->Matrices[i]->message = calloc(dec_ctx->sc->params.size_g, sizeof(GF_ELEMENT*));
        if (dec_ctx->Matrices[i]->message == NULL) {
            fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage\n", fname, i);
            goto AllocError;
        }
        for (int j=0; j<dec_ctx->sc->params.size_g; j++) {
            dec_ctx->Matrices[i]->message[j] = calloc(dec_ctx->sc->params.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->Matrices[i]->message[j] == NULL) {
                fprintf(stderr, "%s: calloc dec_ctx->Matrices[%d]->messsage[%d]\n", fname, i, j);
                goto AllocError;
            }
        }
    }
    if ( (dec_ctx->recent = malloc(sizeof(ID_list))) == NULL ) {
        fprintf(stderr, "%s: malloc dec_ctx->recent", fname);
        goto AllocError;
    }

    dec_ctx->recent->first = dec_ctx->recent->last = NULL;
    memset(dec_ctx->grecent, -1, sizeof(int)*FB_THOLD);             /* set recent decoded generation ids to -1 */
    dec_ctx->newgpos    = 0;
    dec_ctx->grcount    = 0;
    dec_ctx->operations = 0;
    dec_ctx->overhead   = 0;
    dec_ctx->ops1 = dec_ctx->ops2 = 0;
    return dec_ctx;

AllocError:
    free_dec_context_GG(dec_ctx);
    dec_ctx = NULL;
    return NULL;
}

void free_dec_context_GG(struct decoding_context_GG *dec_ctx)
{
    if (dec_ctx == NULL)
        return;
    int i, j, k;
    if (dec_ctx->evolving_checks != NULL) {
        for (i=0; i<dec_ctx->sc->cnum; i++) {
            if (dec_ctx->evolving_checks[i] != NULL)
                free(dec_ctx->evolving_checks[i]);
        }
        free(dec_ctx->evolving_checks);
    }
    if (dec_ctx->check_degrees != NULL)
        free(dec_ctx->check_degrees);
    if (dec_ctx->Matrices != NULL) {
        for (i=0; i<dec_ctx->sc->gnum; i++){
            // Free each decoding matrix
            if (dec_ctx->Matrices[i] != NULL) {
                free_running_matrix(dec_ctx->Matrices[i], dec_ctx->sc->params.size_g);
                dec_ctx->Matrices[i] = NULL;
            }
        }
        free(dec_ctx->Matrices);
    }
    if (dec_ctx->sc != NULL) {
        free_subgen_nbr_list(dec_ctx->sc, gene_nbr);
        snc_free_enc_context(dec_ctx->sc);
    }

    if (dec_ctx->recent != NULL)
        free_list(dec_ctx->recent);
    free(dec_ctx);
    dec_ctx = NULL;
    return;
}

static void free_running_matrix(struct running_matrix *matrix, int rows)
{
    if (matrix != NULL) {
        if (matrix->erased != NULL)
            free(matrix->erased);
        int i;
        if (matrix->row != NULL) {
            for (i=0; i<rows; i++) {
                if (matrix->row[i] != NULL) {
                    if (matrix->row[i]->elem != NULL)
                        free(matrix->row[i]->elem);
                    free(matrix->row[i]);
                }
            }
            free(matrix->row);
        }
        if (matrix->message != NULL) {
            for (i=0; i<rows; i++) {
                if (matrix->message[i] != NULL)
                    free(matrix->message[i]);
            }
            free(matrix->message);
        }
        free(matrix);
        matrix = NULL;
        return;
    }
    return;
}


// cache a new received packet, extract its information and try decode the class it belongs to
void process_packet_GG(struct decoding_context_GG *dec_ctx, struct snc_packet *pkt)
{
    static char fname[] = "process_packet_GG";
    dec_ctx->overhead += 1;

    if (get_loglevel() == TRACE)
        printf("Received: %d g_decoded: %d c_decoded: %d\n", dec_ctx->overhead-1, g_decoded, c_decoded);

    int i, j;
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;

    int gid = pkt->gid;
    if (gid == -1) {
        int pktid = pkt->ucid;
        if (dec_ctx->sc->pp[pktid] != NULL) {
            printf("systematic packet %d has already been decoded/received.\n", pktid);
            return;
        }
        dec_ctx->sc->pp[pktid] = malloc(sizeof(GF_ELEMENT)*pktsize);
        memcpy(dec_ctx->sc->pp[pktid], pkt->syms, sizeof(GF_ELEMENT)*pktsize);
        // Record the as a recently decoded packet
        ID *new_id;
        if ( (new_id = malloc(sizeof(ID))) == NULL )
            fprintf(stderr, "%s: malloc new ID\n", fname);
        new_id->data = pktid;
        new_id->next = NULL;
        append_to_list(dec_ctx->recent, new_id);
        perform_iterative_decoding(dec_ctx);
        return;
    } else {
        // normal GNC packets
        struct running_matrix *matrix = dec_ctx->Matrices[gid];
        if (matrix == NULL)
            return;             // running matrix is freed when the subgeneration is completely decoded

        int DoF_miss = matrix->DoF_miss;
        // 1, Process the packet against the running matrix
        int pivotfound = 0;
        int pivot;
        GF_ELEMENT *pkt_coes = calloc(gensize, sizeof(GF_ELEMENT));
        GF_ELEMENT ce;
        for (i=0; i<gensize; i++) {
            if (dec_ctx->sc->params.gfpower==1) {
                ce = get_bit_in_array(pkt->coes, i);
            } else if (dec_ctx->sc->params.gfpower==8) {
                ce = pkt->coes[i];
            } else {
                ce = read_bits_from_byte_array(pkt->coes, dec_ctx->sc->params.size_g, dec_ctx->sc->params.gfpower, i);
            }
            pkt_coes[i] = ce;
            // if the corresponding source packet has been decoded, remove it from the coded packet
            if (get_bit_in_array(matrix->erased, i) == 1) {
                //find the decoded packet, mask it with this source packet
                int src_id = dec_ctx->sc->gene[gid]->pktid[i];      // index of the corresponding source packet
                galois_multiply_add_region(pkt->syms, dec_ctx->sc->pp[src_id], ce, pktsize);
                pkt_coes[i] = 0;
            }
        }
        // Translate the encoding vector to the sorted form as in the generation
        GF_ELEMENT quotient;
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
            matrix->DoF_miss -= 1;
        }
        free(pkt_coes);
        // 2, Check if the matrix is full rank, if yes, finish decode it
        if (matrix->DoF_miss == 0) {
            decode_generation(dec_ctx, gid); 
            if (get_loglevel() == TRACE)
                printf("Entering perform_iterative_decoding...\n");
            perform_iterative_decoding(dec_ctx);
        }
        if (dec_ctx->finished && get_loglevel() == TRACE) {
            printf("GG splitted operations: %.2f %.2f\n",
                    (double) dec_ctx->ops1/dec_ctx->sc->snum/dec_ctx->sc->params.size_p,
                    (double) dec_ctx->ops2/dec_ctx->sc->snum/dec_ctx->sc->params.size_p);
            printf("Received: %d g_decoded: %d c_decoded: %d\n", dec_ctx->overhead, g_decoded, c_decoded);
        }
        return;
    }
}

// decode packets of a generation via Gaussion elimination
static void decode_generation(struct decoding_context_GG *dec_ctx, int gid)
{
    static char fname[] = "decode_generation";
    //printf("entering decoding_class()...\n");
    struct running_matrix *matrix = dec_ctx->Matrices[gid];
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int i, j, k;
    GF_ELEMENT quotient;
    for (i=gensize-1; i>=0; i--) {
        if (get_bit_in_array(matrix->erased, i) == 1)
            continue;
        /* eliminate all nonzeros above diagonal elements from right to left*/
        for (j=0; j<i; j++) {
            int len = matrix->row[j]->len;
            if (j+len <= i || matrix->row[j]->elem[i-j] == 0)
                continue;
            assert(matrix->row[i]->elem[0]);
            quotient = galois_divide(matrix->row[j]->elem[i-j], matrix->row[i]->elem[0]);
            galois_multiply_add_region(matrix->message[j], matrix->message[i], quotient, pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->ops1 += (pktsize + 1);
            matrix->row[j]->elem[i-j] = 0;
        }
        // transform diagonals to 1
        if (matrix->row[i]->elem[0] != 1) {
            galois_multiply_region(matrix->message[i], galois_divide(1, matrix->row[i]->elem[0]), pktsize);
            dec_ctx->operations += (pktsize + 1);
            dec_ctx->ops1 += (pktsize + 1);
            matrix->row[i]->elem[0] = 1;
        }
        set_bit_in_array(matrix->erased, i);
    }

    // Contruct decoded packets
    int c = 0;                                              // record number of decoded pacekts
    int src_id;
    for (i=0; i<gensize; i++) {
        src_id = dec_ctx->sc->gene[gid]->pktid[i];
        if (dec_ctx->sc->pp[src_id] != NULL) {
            if (get_loglevel() == TRACE)
                printf("%s: packet %d in subgeneration %d is already decoded.\n", fname, src_id, gid);
            continue;
        } else {
            if ( (dec_ctx->sc->pp[src_id] = calloc(pktsize, sizeof(GF_ELEMENT))) == NULL )
                fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, src_id);
            memcpy(dec_ctx->sc->pp[src_id], matrix->message[i], sizeof(GF_ELEMENT)*pktsize);
            if (get_loglevel() == TRACE)
                printf("%s: packet %d decoded from generation %d\n", fname, src_id, gid);
            // Record the decoded packet as a recently decoded packet
            ID *new_id;
            if ( (new_id = malloc(sizeof(ID))) == NULL )
                fprintf(stderr, "%s: malloc new ID\n", fname);
            new_id->data = src_id;
            new_id->next = NULL;
            append_to_list(dec_ctx->recent, new_id);
            c += 1;
        }
    }
    if (get_loglevel() == TRACE)
        printf("%s: %d packets decoded from generation %d\n", fname, c, gid);
    free_running_matrix(matrix, gensize);
    dec_ctx->Matrices[gid] = NULL;

    /* Record recent decoded generations */
    dec_ctx->newgpos = dec_ctx->newgpos % FB_THOLD;
    dec_ctx->grecent[dec_ctx->newgpos] = gid;
    dec_ctx->newgpos++;
    dec_ctx->grcount++;

    g_decoded += c;
}

// This function performs iterative decoding on the precode and GNC code,
// based on the most recently decoded packets from a generation by decode_generation()
static void perform_iterative_decoding(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "perform_iterative_decoding";
    // Perform iterative decoding on the bipartite graph of LDPC precode
    ID *new_decoded = dec_ctx->recent->first;
    while (new_decoded != NULL) {
        int new_id = new_decoded->data;
        if (new_id >= dec_ctx->sc->snum) {
            new_decoded_check_packet(dec_ctx, new_id);
        } else {
            new_decoded_source_packet(dec_ctx, new_id);
            if (dec_ctx->finished)
                return;     // Leave iterative decoding if all source packets are decoded
        }
        check_for_new_recoverables(dec_ctx);
        new_decoded = new_decoded->next;
    }

    // Perform iterative generation decoding
    update_generations(dec_ctx);
    int new_decodable_gid = check_for_new_decodables(dec_ctx);
    if (new_decodable_gid != -1) {
        decode_generation(dec_ctx, new_decodable_gid);
        perform_iterative_decoding(dec_ctx);
    }
}

// Precedures to take when a source packet is decoded from a generation
static void new_decoded_source_packet(struct decoding_context_GG *dec_ctx, int pkt_id)
{
    static char fname[] = "new_decoded_source_packet";
    dec_ctx->decoded   += 1;
    dec_ctx->originals += 1;
    if (dec_ctx->originals == dec_ctx->sc->snum) {
        dec_ctx->finished = 1;
        return;      // Don't need to worry about the remaining check packets
    }
    // If the decoded packet is an original source packet back substitute it 
    // into those known LDPC check packets
    if (dec_ctx->sc->graph == NULL) {
        //no precode
        return;
    }
    NBR_node *nb = dec_ctx->sc->graph->r_nbrs_of_l[pkt_id]->first;
    while (nb != NULL) {
        int check_id = nb->data;
        // If the corresponding check packet is not yet decoded, the evolving packet area
        // and the corresponding degree can be used to record the evolution of the packet.
        if (dec_ctx->evolving_checks[check_id] == NULL) {
            dec_ctx->evolving_checks[check_id] = calloc(dec_ctx->sc->params.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->evolving_checks[check_id] == NULL)
                fprintf(stderr, "%s: calloc evolving_checks[%d]\n", fname, check_id);
        }
        // mask information bits
        galois_multiply_add_region(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], nb->ce, dec_ctx->sc->params.size_p);
        dec_ctx->operations += dec_ctx->sc->params.size_p;
        dec_ctx->ops2 += dec_ctx->sc->params.size_p;
        dec_ctx->check_degrees[check_id] -= 1;    // reduce check packet's unknown degree
        // We don't actually need to remove nodes from bipartite during iterative decoding. We already track
        // the degree changes in evolving_packets and check_degrees. Leaving sc->graph intact is a better choice.
        // if (remove_from_list(dec_ctx->sc->graph->l_nbrs_of_r[check_id], pkt_id) == -1)
        //    fprintf(stderr, "%s: remove %d from l_nbrs_of_r[%d]\n", fname, pkt_id, check_id);

        nb = nb->next;
    }
}

// Procedures to take when a check packet is decoded from a generation
static void new_decoded_check_packet(struct decoding_context_GG *dec_ctx, int pkt_id)
{
    static char fname[] = "new_decoded_check_packet";
    dec_ctx->decoded += 1;

    int check_id = pkt_id - dec_ctx->sc->snum;                         // ID of check packet
    if (get_loglevel() == TRACE)
        printf("%s: degree of new decoded check %d is %d\n", \
                fname, pkt_id, dec_ctx->check_degrees[check_id]);
    if (dec_ctx->evolving_checks[check_id] == NULL) {
        // Evolving area is empty, meaning that no source neighbors of the check packet
        // has been decoded yet. So make a copy of the decoded check packet for later evolving
        dec_ctx->evolving_checks[check_id] = calloc(dec_ctx->sc->params.size_p, sizeof(GF_ELEMENT));
        if (dec_ctx->evolving_checks[check_id] == NULL)
            fprintf(stderr, "%s: calloc evolving_checks[%d]\n", fname, check_id);
        memcpy(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], sizeof(GF_ELEMENT)*dec_ctx->sc->params.size_p);
    } else {
        // Some source neighbors have been decoded and therefore updated the evolving area,
        // we need to mask the actual check packet content against the evolving area.
        galois_multiply_add_region(dec_ctx->evolving_checks[check_id], dec_ctx->sc->pp[pkt_id], 1, dec_ctx->sc->params.size_p);
        dec_ctx->operations += dec_ctx->sc->params.size_p;
        dec_ctx->ops2 += dec_ctx->sc->params.size_p;
    }
}

// This function is part of iterative LDPC _precode_ decoding, which
// checks for new recoverable source/check packet from parity-check packets.
static int check_for_new_recoverables(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "check_for_new_recoverables";
    int snum = dec_ctx->sc->snum;
    int has_new_recoverable = -1;
    // check each check node
    for (int i=0; i<dec_ctx->sc->cnum; i++) {
        if (dec_ctx->check_degrees[i] == 1
                && dec_ctx->sc->pp[i+snum] != NULL
                && !exist_in_list(dec_ctx->recent, i+snum)) {
            // The check packet is already decoded from some previous generations and its degree is
            // reduced to 1, meaning that it connects to a unrecovered source neighboer. Let's find and recover this
            // source neighbor. Note that there is at most 1 (i.e., 0 or 1) packet may be recovered from each check.
            // So it is safe to break whenever a NULL pp is found (i.e., a packet not recovered yet)
            int found = 0;
            int pktid;
            NBR_node *nb = dec_ctx->sc->graph->l_nbrs_of_r[i]->first;
            while (nb != NULL) {
                pktid = nb->data;
                if (dec_ctx->sc->pp[pktid] == NULL) {
                    found = 1;
                    break;
                }
                nb = nb->next;
            }
            if (found) {
                GF_ELEMENT ce = nb->ce;
                if (get_loglevel() == TRACE)
                    printf("%s: source packet %d is recovered from check %d\n", fname, pktid, i+snum);
                dec_ctx->sc->pp[pktid] = calloc(dec_ctx->sc->params.size_p, sizeof(GF_ELEMENT));
                if (dec_ctx->sc->pp[pktid] == NULL)
                    fprintf(stderr, "%s: calloc sc->pp[%d]\n", fname, pktid);
                if (ce == 1)
                    memcpy(dec_ctx->sc->pp[pktid], dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT)*dec_ctx->sc->params.size_p);
                else {
                    galois_multiply_add_region(dec_ctx->sc->pp[pktid], dec_ctx->evolving_checks[i], galois_divide(1, ce), dec_ctx->sc->params.size_p);
                    dec_ctx->operations += dec_ctx->sc->params.size_p + 1;
                    dec_ctx->ops2 += dec_ctx->sc->params.size_p + 1;
                }
                // Record the decoded packet as a recently decoded packet
                ID *new_id;
                if ( (new_id = malloc(sizeof(ID))) == NULL )
                    fprintf(stderr, "%s: malloc new ID\n", fname);
                new_id->data = pktid;
                new_id->next = NULL;
                append_to_list(dec_ctx->recent, new_id);

                c_decoded += 1;
            }
            dec_ctx->check_degrees[i] = 0;
        }
        if (dec_ctx->sc->pp[i+snum] == NULL && dec_ctx->check_degrees[i] == 0) {
            // The check packet is recovered because all of its source neighbors are known, so it is recoverable
            if (get_loglevel() == TRACE)
                printf("%s: check packet %d is recoverable after all the connecting source neighbors are known\n", fname, i+snum);
            dec_ctx->sc->pp[i+snum] = calloc(dec_ctx->sc->params.size_p, sizeof(GF_ELEMENT));
            if (dec_ctx->sc->pp[i+snum] == NULL)
                fprintf(stderr, "%s: calloc sc->pp[%d]", fname, i+snum);
            memcpy(dec_ctx->sc->pp[i+snum], dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT)*dec_ctx->sc->params.size_p);
            // Record a recently decoded packet
            ID *new_id;
            if ( (new_id = malloc(sizeof(ID))) == NULL )
                fprintf(stderr, "%s: malloc new ID\n", fname);

            new_id->data = i+snum;
            new_id->next = NULL;
            append_to_list(dec_ctx->recent, new_id);

            c_decoded += 1;
        }

    }
    return has_new_recoverable;
}

// Update non-decoded generations with recently decoded packets
static void update_generations(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "update_generations";

    ID *precent = dec_ctx->recent->first;  // pointer to index of recently decoded packets
    while (precent != NULL) {
        int src_id = precent->data;
        // Check all generations that contain this source packet
        ID *item = gene_nbr[src_id]->first;
        while (item != NULL) {
            int gid = item->data;
            item = item->next;
            if (dec_ctx->Matrices[gid] == NULL || dec_ctx->Matrices[gid]->DoF_miss == 0)
                continue;
            int pos = has_item(dec_ctx->sc->gene[gid]->pktid, src_id, dec_ctx->sc->params.size_g);
            if (pos != -1
                && get_bit_in_array(dec_ctx->Matrices[gid]->erased, pos) == 0) {
                // The recently decoded packet is also in this generation
                long ops = update_running_matrix(dec_ctx, gid, src_id, pos);
                dec_ctx->operations += ops;
                dec_ctx->ops2 += ops;
            }
        }
        precent = precent->next;
    }
    // Clean up recently decoded packet ID list
    clear_list(dec_ctx->recent);
}


// Update running matrix with a decoded source packet
static long update_running_matrix(struct decoding_context_GG *dec_ctx, int gid, int sid, int index)
{
    static char fname[] = "update_running_matrix";
    //printf("entering update_running_matrix()...\n");
    
    // For each running matrix containing the just decoded packet [sid], a singleton vector will be
    // placed in the matrix at the index's row, which corresponds to the position of the packet in 
    // this subgeneration. All the nonzero vectors above this row will remove the corresponding column
    // (i.e., set coding coefficient to 0). If there was a nonzero vector at the index's row, it will
    // be replaced, and the previous row will be processed against rows below it, until it is placed
    // in some zero vectors or being reduced to zero and discarded.
    int i, j, k;
    long operations = 0;
    struct running_matrix *matrix = dec_ctx->Matrices[gid];
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;

    // 1) Process vectors above the index's row
    GF_ELEMENT ce;
    for (i=0; i<index; i++) {
        if (matrix->row[i] != NULL 
            && matrix->row[i]->len > (index-i)
            && matrix->row[i]->elem[index-i] !=0) {
            ce = matrix->row[i]->elem[index-i];
            galois_multiply_add_region(matrix->message[i], dec_ctx->sc->pp[sid], ce, pktsize);
            matrix->row[i]->elem[index-i] = 0;
        }
    }
    
    // 2) Place singleton vector to the index's row
    GF_ELEMENT *coes = calloc(gensize, sizeof(GF_ELEMENT));
    GF_ELEMENT *syms = calloc(pktsize, sizeof(GF_ELEMENT));
    // Process the packet just as a normal received packet of the matrix
    if (matrix->row[index] == NULL) {
        // Just fill in the row with a singleton row and the decoded packet
        matrix->row[index] = (struct row_vector*) malloc(sizeof(struct row_vector));
        if (matrix->row[index] == NULL)
            fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, index);
        matrix->row[index]->len = 1;
        matrix->row[index]->elem = (GF_ELEMENT *) calloc(matrix->row[index]->len, sizeof(GF_ELEMENT));
        if (matrix->row[index]->elem == NULL)
            fprintf(stderr, "%s: calloc matrix->row[%d]->elem failed\n", fname, index);
        matrix->row[index]->elem[0] = 1;
        memcpy(matrix->message[index], dec_ctx->sc->pp[sid],  pktsize*sizeof(GF_ELEMENT));
        matrix->DoF_miss -= 1;
    } else {
        // Let's do swap first, and leave a singleton vector at the row 
        memcpy(coes+index, matrix->row[index]->elem, sizeof(GF_ELEMENT)*(gensize-index));
        memset(matrix->row[index]->elem, 0, sizeof(GF_ELEMENT)*(gensize-index));
        matrix->row[index]->len = 1;
        matrix->row[index]->elem[0] = 1;
        memcpy(syms, matrix->message[index], sizeof(GF_ELEMENT)*pktsize);
        memcpy(matrix->message[index], dec_ctx->sc->pp[sid], sizeof(GF_ELEMENT)*pktsize);
        // and then process the previous row as if it was a received vector
        // process it against rows below the current row
        int pivotfound = 0;
        int pivot;
        GF_ELEMENT quotient;
        for (i=index; i<gensize; i++) {
            if (coes[i] != 0) {
                if (matrix->row[i] != NULL) {
                    quotient = galois_divide(coes[i], matrix->row[i]->elem[0]);
                    galois_multiply_add_region(&(coes[i]), matrix->row[i]->elem, quotient, matrix->row[i]->len);
                    galois_multiply_add_region(syms, matrix->message[i], quotient, pktsize);
                    dec_ctx->operations += 1 + matrix->row[i]->len + pktsize;
                    dec_ctx->ops1 += 1 + matrix->row[i]->len + pktsize;
                } else {
                    pivotfound = 1;
                    pivot = i;
                    break;
                }
            }
        }
        // if pivot found, save it
        if (pivotfound) {
            matrix->row[pivot] = (struct row_vector*) malloc(sizeof(struct row_vector));
            if (matrix->row[pivot] == NULL)
                fprintf(stderr, "%s: malloc dec_ctx->row[%d] failed\n", fname, pivot);
            matrix->row[pivot]->len = gensize - pivot;
            matrix->row[pivot]->elem = (GF_ELEMENT *) calloc(matrix->row[pivot]->len, sizeof(GF_ELEMENT));
            if (matrix->row[pivot]->elem == NULL)
                fprintf(stderr, "%s: calloc matrix->row[%d]->elem failed\n", fname, pivot);
            memcpy(matrix->row[pivot]->elem, &(coes[pivot]), matrix->row[pivot]->len*sizeof(GF_ELEMENT));
            memcpy(matrix->message[pivot], syms,  pktsize*sizeof(GF_ELEMENT));
            matrix->DoF_miss -= 1;
        }
    }
    free(coes);
    free(syms);
    set_bit_in_array(matrix->erased, index);       // mark the corresponding column as erased
    return operations;
}

// Check if there is new decodable generations
// Return:
//    Decodable generation id, or else -1.
static int check_for_new_decodables(struct decoding_context_GG *dec_ctx)
{
    static char fname[] = "check_for_new_decodables";
    int i, j, k;
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        struct running_matrix *matrix = dec_ctx->Matrices[i];
        if (matrix != NULL && matrix->DoF_miss == 0)
            return i;
    }
    return -1;
}


/**
 * Save a decoding context to a file
 * Return values:
 *   On success: bytes written
 *   On error: -1
 */
long save_dec_context_GG(struct decoding_context_GG *dec_ctx, const char *filepath)
{
    long filesize = 0;
    int d_type = GG_DECODER;
    FILE *fp;
    if ((fp = fopen(filepath, "w")) == NULL) {
        fprintf(stderr, "Cannot open %s to save decoding context\n", filepath);
        return (-1);
    }
    int i, j, k;
    int gensize = dec_ctx->sc->params.size_g;
    int pktsize = dec_ctx->sc->params.size_p;
    int numpp   = dec_ctx->sc->snum + dec_ctx->sc->cnum;
    // Write snc params
    filesize += fwrite(&dec_ctx->sc->params, sizeof(struct snc_parameters), 1, fp);
    // Write decoder type
    filesize += fwrite(&(d_type), sizeof(int), 1, fp);
    // Save already decoded packets in dec_ctx->sc->pp
    filesize += fwrite(&dec_ctx->decoded, sizeof(int), 1, fp);
    for (i=0; i<numpp; i++) {
        if (dec_ctx->sc->pp[i] != NULL) {
            filesize += fwrite(&i, sizeof(int), 1, fp);  // pktid
            filesize += fwrite(dec_ctx->sc->pp[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    // Save evolving check packets
    int count = 0;
    for (i=0; i<dec_ctx->sc->cnum; i++) {
        if (dec_ctx->evolving_checks[i] != NULL)
            count++;
    }
    filesize += fwrite(&count, sizeof(int), 1, fp);  // evolving packets non-NULL count
    for (i=0; i<dec_ctx->sc->cnum; i++) {
        if (dec_ctx->evolving_checks[i] != NULL) {
            filesize += fwrite(&i, sizeof(int), 1, fp);  // check id
            filesize += fwrite(dec_ctx->evolving_checks[i], sizeof(GF_ELEMENT), pktsize, fp);
        }
    }
    // Save check degrees
    filesize += fwrite(dec_ctx->check_degrees, sizeof(int), dec_ctx->sc->cnum, fp);
    filesize += fwrite(&dec_ctx->finished, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->decoded, sizeof(int), 1, fp);  // Yes, I know. I saved this value twice!
    filesize += fwrite(&dec_ctx->originals, sizeof(int), 1, fp);
    // Save running matrices
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        // Use a flag to indicate whether the matrix is freed
        int mflag = dec_ctx->Matrices[i] == NULL ? 0 : 1;
        filesize += fwrite(&mflag, sizeof(int), 1, fp);
        if (mflag == 0)
            continue;    // pass NULL matrices

        filesize += fwrite(&dec_ctx->Matrices[i]->DoF_miss, sizeof(int), 1, fp);
        int nflags = ALIGN(dec_ctx->sc->params.size_g, 8);
        filesize += fwrite(dec_ctx->Matrices[i]->erased, 1, nflags, fp);
        // Write rows of running matrix
        for (j=0; j<gensize; j++) {
            // Use a flag to indicate whether a row exist
            int rowlen = dec_ctx->Matrices[i]->row[j] == NULL ? 0 : dec_ctx->Matrices[i]->row[j]->len;
            filesize += fwrite(&rowlen, sizeof(int), 1, fp);
            if (rowlen != 0) {
                filesize += fwrite(dec_ctx->Matrices[i]->row[j]->elem, sizeof(GF_ELEMENT), rowlen, fp);
                filesize += fwrite(dec_ctx->Matrices[i]->message[j], sizeof(GF_ELEMENT), pktsize, fp);
            }
        }
    }
    // Save recent ID_list
    count = 0;
    ID *id = dec_ctx->recent->first;
    while (id != NULL) {
        count++;
        id = id->next;
    }
    filesize += fwrite(&count, sizeof(int), 1, fp);  // Number of recent IDs in the list
    id = dec_ctx->recent->first;
    while (count > 0) {
        fwrite(&id->data, sizeof(int), 1, fp);
        count--;
        id = id->next;
    }
    // Save performance index
    filesize += fwrite(&dec_ctx->overhead, sizeof(int), 1, fp);
    filesize += fwrite(&dec_ctx->operations, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops1, sizeof(long long), 1, fp);
    filesize += fwrite(&dec_ctx->ops2, sizeof(long long), 1, fp);
    fclose(fp);
    return filesize;
}

struct decoding_context_GG *restore_dec_context_GG(const char *filepath)
{
    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL) {
        fprintf(stderr, "Cannot open %s to load decoding context\n", filepath);
        return NULL;
    }
    struct snc_parameters sp;
    fread(&sp, sizeof(struct snc_parameters), 1, fp);
    // Create a fresh decoding context
    struct decoding_context_GG *dec_ctx = create_dec_context_GG(&sp);
    if (dec_ctx == NULL) {
        fprintf(stderr, "malloc decoding_context_GG failed\n");
        return NULL;
    }
    
    build_subgen_nbr_list(dec_ctx->sc);

    // Restore decoding context from file
    fseek(fp, sizeof(int), SEEK_CUR);  // skip decoding_type field
    // Restore already decoded packets
    int i, j, k;
    fread(&dec_ctx->decoded, sizeof(int), 1, fp);
    for (i=0; i<dec_ctx->decoded; i++) {
        int pktid;
        fread(&pktid, sizeof(int), 1, fp);
        dec_ctx->sc->pp[pktid] = calloc(sp.size_p, sizeof(GF_ELEMENT));
        fread(dec_ctx->sc->pp[pktid], sizeof(GF_ELEMENT), sp.size_p, fp);
    }
    // Restore evolving packets
    int count;
    fread(&count, sizeof(int), 1, fp);
    for (i=0; i<count; i++) {
        int evoid;
        fread(&evoid, sizeof(int), 1, fp);
        dec_ctx->evolving_checks[evoid] = calloc(sp.size_p, sizeof(GF_ELEMENT));
        fread(dec_ctx->evolving_checks[evoid], sizeof(GF_ELEMENT), sp.size_p, fp);
    }
    // Restore check degrees
    fread(dec_ctx->check_degrees, sizeof(int), dec_ctx->sc->cnum, fp);
    fread(&dec_ctx->finished, sizeof(int), 1, fp);
    fread(&dec_ctx->decoded, sizeof(int), 1, fp);
    fread(&dec_ctx->originals, sizeof(int), 1, fp);
    // Restore running matrices
    // Note that running matrices' memory were already allocated in creating_dec_context
    // If a matrix was freed in the previous stored decoder context, we need to free again 
    // after restoring it.
    int mflag;  // flag indicating whether a matrix was already freed when saving the decoder context
    for (i=0; i<dec_ctx->sc->gnum; i++) {
        fread(&mflag, sizeof(int), 1, fp);
        if (mflag == 0) {
            free_running_matrix(dec_ctx->Matrices[i], sp.size_g);
            continue;
        }
        fread(&dec_ctx->Matrices[i]->DoF_miss, sizeof(int), 1, fp);
        int nflags = ALIGN(dec_ctx->sc->params.size_g, 8);
        fread(dec_ctx->Matrices[i]->erased, 1, nflags, fp);
        // read rows
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
    // Restore recent ID_list
    fread(&count, sizeof(int), 1, fp);
    ID *new_id;
    for (i=0; i<count; i++) {
        if ( (new_id = malloc(sizeof(ID))) == NULL )
            fprintf(stderr, "malloc new ID failed\n");
        fread(&new_id->data, sizeof(int), 1, fp);
        new_id->next = NULL;
        append_to_list(dec_ctx->recent, new_id);
    }
    // Restore performance index
    fread(&dec_ctx->overhead, sizeof(int), 1, fp);
    fread(&dec_ctx->operations, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops1, sizeof(long long), 1, fp);
    fread(&dec_ctx->ops2, sizeof(long long), 1, fp);
    fclose(fp);
    return dec_ctx;
}


