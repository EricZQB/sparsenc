/**************************************************************
 * sncEncoder.c
 *
 * Functions for SNC encoding. Coded packets can be generated
 * from memory buffer or files.
 **************************************************************/
#include <math.h>
#include <sys/time.h>
#include "common.h"
#include "galois.h"
#include "sparsenc.h"

static int create_context_from_params(struct snc_context *sc);
static int verify_code_parameter(struct snc_parameters *sp);
static void perform_precoding(struct snc_context *sc);
static int group_packets_rand(struct snc_context *sc);
static int group_packets_pseudorand(struct snc_context *sc);
static int group_packets_band(struct snc_context *sc);
static int group_packets_windwrap(struct snc_context *sc);
static void encode_packet(struct snc_context *sc, int gid, struct snc_packet *pkt);
static int schedule_generation(struct snc_context *sc);
static int banded_nonuniform_sched(struct snc_context *sc);
/*
 * Create a GNC context containing meta information about the data to be encoded.
 *   buf      - Buffer containing bytes of data to be encoded
 *   sc       - pointer to the snc_context where the context will be stored
 *   snc_context is defined in src/common.h
 *
 * Return
 *   0  - Create successfully
 *   -1 - Create failed
 */
struct snc_context *snc_create_enc_context(unsigned char *buf, struct snc_parameters *sp)
{
    static char fname[] = "snc_create_enc_context";
    // Set log level
    char *loglevel = getenv("SNC_LOG_LEVEL");
    if (loglevel != NULL)
        set_loglevel(loglevel);

    // Allocate file_context
    struct snc_context *sc;
    if ( (sc = calloc(1, sizeof(struct snc_context))) == NULL ) {
        fprintf(stderr, "%s: calloc file_context\n", fname);
        return NULL;
    }
    sc->params.datasize = sp->datasize;
    sc->params.size_p   = sp->size_p;
    sc->params.size_c   = sp->size_c;
    sc->params.size_b   = sp->size_b;
    sc->params.size_g   = sp->size_g;
    sc->params.type     = sp->type;
    sc->params.bpc      = sp->bpc;
    sc->params.bnc      = sp->bnc;
    sc->params.sys      = sp->sys;
    sc->params.seed     = sp->seed;
    /* Seed local random number generator for precoding and/or random grouping
     *
     *   If creating a completely new snc context, seed is -1 by default. We
     *   will seed using current time stamp.
     *
     *   In other cases, we might create a sc based on sp that is in a file or
     *   received from a network. In this case, seed != -1 and we just use the
     *   given seed.
     */
    if (sc->params.seed == -1) {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        sc->params.seed = tv.tv_sec * 1000 + tv.tv_usec / 1000; // seed use microsec
    }
    init_genrand(sc->params.seed);
    sp->seed = sc->params.seed;  // set seed in the passed-in argument as well
    // Determine packet and generation numbers
    int num_src = ALIGN(sc->params.datasize, sc->params.size_p);
    int num_chk = sp->size_c;
    sc->snum  = num_src;  // Number of source packets
    sc->cnum  = num_chk;  // Number of check packets
    if (sc->params.type == BAND_SNC) {
        sc->gnum  = ALIGN((num_src+num_chk-sc->params.size_g), sc->params.size_b) + 1;
    } else {
        sc->gnum  = ALIGN( (num_src+num_chk), sc->params.size_b);
    }
    /*
     * Verify code parameter
     */
    if (verify_code_parameter(&sc->params) != 0) {
        fprintf(stderr, "%s: code parameter is invalid.\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }
    /*
     * Create generations, bipartite graph
     */
    if (create_context_from_params(sc) != 0) {
        fprintf(stderr, "%s: create_context_from_params\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }

    // Allocating pointers to data
    if ((sc->pp = calloc(sc->snum+sc->cnum, sizeof(GF_ELEMENT*))) == NULL) {
        fprintf(stderr, "%s: calloc sc->pp\n", fname);
        snc_free_enc_context(sc);
        return NULL;
    }

    constructField();   // Construct Galois Field
    if (buf != NULL) {
        int alread = 0;
        int i;
        // Load source packets
        for (i=0; i<sc->snum; i++) {
            sc->pp[i] = calloc(sc->params.size_p, sizeof(GF_ELEMENT));
            int toread = (alread+sc->params.size_p) <= sc->params.datasize ? sc->params.size_p : sc->params.datasize-alread;
            memcpy(sc->pp[i], buf+alread, toread*sizeof(GF_ELEMENT));
            alread += toread;
        }
        // Allocate parity-check packet space
        for (i=0; i<sc->cnum; i++)
            sc->pp[sc->snum+i] = calloc(sc->params.size_p, sizeof(GF_ELEMENT));
        perform_precoding(sc);
    }

    return sc;
}

inline struct snc_parameters *snc_get_parameters(struct snc_context *sc)
{
    if (sc == NULL) {
        return NULL;
    } else {
        return &(sc->params);
    }
}

/*
 * Load data from file  into a snc_context.
 *   start - starting point to read file
 * It is caller's responsibility to ensure:
 *   start + sc->params.datasize <= fileSize
 */
int snc_load_file_to_context(const char *filepath, long start, struct snc_context *sc)
{
    static char fname[] = "snc_load_file_to_context";

    FILE *fp;
    if ((fp = fopen(filepath, "r")) == NULL)
        return (-1);
    fseek(fp, 0, SEEK_END);
    if ((ftell(fp) - start) < sc->params.datasize) {
        return (-1);
    }
    fseek(fp, start, SEEK_SET);  // seek to position start
    int alread = 0;
    int i;
    for (i=0; i<sc->snum; i++) {
        sc->pp[i] = calloc(sc->params.size_p, sizeof(GF_ELEMENT));
        int toread = (alread+sc->params.size_p) <= sc->params.datasize ? sc->params.size_p : sc->params.datasize-alread;
        if (fread(sc->pp[i], sizeof(GF_ELEMENT), toread, fp) != toread) {
            fprintf(stderr, "%s: fread sc->pp[%d]\n", fname, i);
            return (-1);
        }
        alread += toread;
    }
    fclose(fp);
    // Allocate parity-check packet space
    for (i=0; i<sc->cnum; i++)
        sc->pp[sc->snum+i] = calloc(sc->params.size_p, sizeof(GF_ELEMENT));
    perform_precoding(sc);
    return (0);
}

static int verify_code_parameter(struct snc_parameters *sp)
{
    if (sp->size_b > sp->size_g) {
        fprintf(stderr, "code spmeter error: size_b > size_g\n");
        return(-1);
    }
    if (sp->size_g*sp->size_p > sp->datasize) {
        fprintf(stderr, "code spmeter error: size_b X size_p > datasize\n");
        return(-1);
    }
    return(0);
}

/*
 * Create snc context using parameters
 */
static int create_context_from_params(struct snc_context *sc)
{
    static char fname[] = "snc_create_enc_context_params";
    // Inintialize generation structures
    sc->gene  = calloc(sc->gnum, sizeof(struct subgeneration*));
    if ( sc->gene == NULL ) {
        fprintf(stderr, "%s: malloc sc->gene\n", fname);
        return(-1);
    }
    for (int j=0; j<sc->gnum; j++) {
        sc->gene[j] = malloc(sizeof(struct subgeneration));
        if ( sc->gene[j] == NULL ) {
            fprintf(stderr, "%s: malloc sc->gene[%d]\n", fname, j);
            return(-1);
        }
        sc->gene[j]->gid = -1;
        sc->gene[j]->pktid = malloc(sizeof(int)*sc->params.size_g);       // Use malloc because pktid needs to be initialized as -1's later
        if ( sc->gene[j]->pktid == NULL ) {
            fprintf(stderr, "%s: malloc sc->gene[%d]->pktid\n", fname, j);
            return(-1);
        }
        memset(sc->gene[j]->pktid, -1, sizeof(int)*sc->params.size_g);
    }
    sc->nccount = calloc(sc->gnum, sizeof(int));
    if (sc->nccount == NULL) {
        fprintf(stderr, "%s: calloc sc->nccount\n", fname);
        return (-1);
    }
    sc->count = 0;

    int coverage;
    if (sc->params.type == RAND_SNC) {
        coverage = group_packets_rand(sc);
        //coverage = group_packets_pseudorand(sc);
    } else if (sc->params.type == BAND_SNC) {
        coverage = group_packets_band(sc);
    } else if (sc->params.type == WINDWRAP_SNC) {
        coverage = group_packets_windwrap(sc);
    }
    // Creating bipartite graph of the precode
    if (sc->cnum != 0) {
        if ( (sc->graph = malloc(sizeof(BP_graph))) == NULL ) {
            fprintf(stderr, "%s: malloc BP_graph\n", fname);
            return (-1);
        }
        sc->graph->binaryce = sc->params.bpc;     // If precode in GF(2), edges use 1 as coefficient
        if (create_bipartite_graph(sc->graph, sc->snum, sc->cnum) < 0)
            return (-1);
    }
    return(0);
}

void snc_free_enc_context(struct snc_context *sc)
{
    if (sc == NULL)
        return;
    int i;
    if (sc->pp != NULL) {
        for (i=sc->snum+sc->cnum-1; i>=0; i--) {
            if (sc->pp[i] != NULL) {
                free(sc->pp[i]);
                sc->pp[i] = NULL;
            }
        }
        free(sc->pp);
    }
    if (sc->gene != NULL) {
        for (i=sc->gnum-1; i>=0; i--) {
            free(sc->gene[i]->pktid);  // free packet IDs
            free(sc->gene[i]);         // free generation itself
            sc->gene[i] = NULL;
        }
        free(sc->gene);
    }
    if (sc->graph != NULL)
        free_bipartite_graph(sc->graph);
    if (sc->nccount != NULL)
        free(sc->nccount);
    free(sc);
    sc = NULL;
    return;
}

unsigned char *snc_recover_data(struct snc_context *sc)
{
    static char fname[] = "snc_recover_data";
    long datasize = sc->params.datasize;
    long alwrote = 0;
    long towrite = datasize;

    unsigned char *data;
    if ( (data = malloc(datasize)) == NULL) {
        fprintf(stderr, "%s: malloc(datasize) failed.\n", fname);
        return NULL;
    }
    memset(data, 0, sizeof(unsigned char)*datasize);
    int pc = 0;
    while (alwrote < datasize) {
        towrite = ((alwrote + sc->params.size_p) <= datasize) ? sc->params.size_p : datasize - alwrote;
        memcpy(data+alwrote, sc->pp[pc++], sizeof(GF_ELEMENT)*towrite);
        alwrote += towrite;
    }
    return data;
}

// Wrapper of free()
void snc_free_recovered(unsigned char *data)
{
    if (data != NULL)
        free(data);
    return;
}

/**
 * Recover data to file.
 * fp has to be opened in 'a' (append) mode
 **/
long snc_recover_to_file(const char *filepath, struct snc_context *sc)
{
    static char fname[] = "snc_recover_to_file";
    long datasize = sc->params.datasize;
    long alwrote = 0;
    long towrite = datasize;

    if (get_loglevel() == TRACE)
        printf("Writing to decoded file.\n");
    FILE *fp;
    if ((fp = fopen(filepath, "a")) == NULL)
        return (-1);

    int pc = 0;
    while (alwrote < datasize) {
        towrite = ((alwrote + sc->params.size_p) <= datasize) ? sc->params.size_p : datasize - alwrote;
        if (fwrite(sc->pp[pc], sizeof(GF_ELEMENT), towrite, fp) != towrite)
            fprintf(stderr, "%s: fwrite sc->pp[%d]\n", fname, pc);
        pc++;
        alwrote += towrite;
    }
    fclose(fp);
    return alwrote;
}

// perform systematic LDPC precoding against SRC pkt list and results in a LDPC pkt list
static void perform_precoding(struct snc_context *sc)
{
    static char fname[] = "perform_precoding";

    int i, j;
    for (i=0; i<sc->cnum; i++) {
        // Encoding check packet according to the LDPC graph
        NBR_node *nb = sc->graph->l_nbrs_of_r[i]->first;
        while(nb != NULL) {
            int sid = nb->data;  // index of source packet
            // XOR information content
            galois_multiply_add_region(sc->pp[i+sc->snum], sc->pp[sid], nb->ce, sc->params.size_p);
            // move to next possible neighbour node of current check
            nb = nb->next;
        }
    }
}

/*
 * This routine uses a deterministic grouping scheme, so the need of sending
 * grouping information to clients is removed. The only information clients
 * need to know is the number of packets, base size, and generation size.
 */
static int group_packets_pseudorand(struct snc_context *sc)
{
    int num_p = sc->snum + sc->cnum;
    int num_g = sc->gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int rotate = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        // split packets into disjoint groups
        for (j=0; j<sc->params.size_b; j++) {
            index = (i * sc->params.size_b + j) % num_p;  // source packet index

            while (has_item(sc->gene[i]->pktid, index, j) != -1)
                index++;
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }

        // fill in the rest of the generation with packets from other generations
        for (j=sc->params.size_b; j<sc->params.size_g; j++) {
            int magicX = sc->params.size_b + num_g - sc->params.size_g;
            magicX = 7 <= magicX ? 7 : magicX;      // Make sure start won't be negative
            int tmp = i - (j - sc->params.size_b + magicX);
            int start = tmp >= 0 ? tmp : tmp+num_g;
            if (start == i)
                start++;
            index = (start * sc->params.size_b + (j - sc->params.size_b + rotate) % (sc->params.size_g)) % num_p;
            while (has_item(sc->gene[i]->pktid, index, j) != -1) {
                index++;
                index = index % num_p;
            }
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }
        rotate = (rotate + 7) % (sc->params.size_g);
    }
    int coverage = 0;
    for (i=0; i<num_p; i++)
        coverage += selected[i];

    free(selected);
    return coverage;
}


/*
 * Use local RNG to group packets
 */
static int group_packets_rand(struct snc_context *sc)
{
    int num_p = sc->snum + sc->cnum;
    int num_g = sc->gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int rotate = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        // split packets into disjoint groups
        for (j=0; j<sc->params.size_b; j++) {
            index = (i * sc->params.size_b + j) % num_p;  // source packet index

            while (has_item(sc->gene[i]->pktid, index, j) != -1)
                index = genrand_int32() % num_p;
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }

        // fill in the rest of the generation with packets from other generations
        for (j=sc->params.size_b; j<sc->params.size_g; j++) {
            index = genrand_int32() % num_p;
            while (has_item(sc->gene[i]->pktid, index, j) != -1) {
                index = genrand_int32() % num_p;
            }
            sc->gene[i]->pktid[j] = index;
            selected[index] += 1;
        }
    }
    int coverage = 0;
    for (i=0; i<num_p; i++)
        coverage += selected[i];

    free(selected);
    return coverage;
}

/*
 * Group packets to generations that overlap head-to-toe. Each generation's
 * encoding coefficients form a band in GDM.
 */
static int group_packets_band(struct snc_context *sc)
{
    int num_p = sc->snum + sc->cnum;
    int num_g = sc->gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int leading_pivot = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        leading_pivot = i * sc->params.size_b;
        if (leading_pivot > num_p - sc->params.size_g) {
            if (get_loglevel() == TRACE)
                printf("Band lead of gid: %d is modified\n", i);
            leading_pivot = num_p - sc->params.size_g;
        }
        for (j=0; j<sc->params.size_g; j++) {
            index = leading_pivot + j;
            selected[index] += 1;
            sc->gene[i]->pktid[j] = index;
        }
    }
    int coverage = 0;
    for (i=0; i<num_p; i++)
        coverage += selected[i];

    free(selected);
    return coverage;
}

/*
 * Group packets to generations that overlap consecutively. Wrap around if needed.
 */
static int group_packets_windwrap(struct snc_context *sc)
{
    int num_p = sc->snum + sc->cnum;
    int num_g = sc->gnum;

    int *selected = calloc(num_p, sizeof(int));

    int i, j;
    int index;
    int leading_pivot = 0;
    for (i=0; i<num_g; i++) {
        sc->gene[i]->gid = i;
        leading_pivot = i * sc->params.size_b;
        for (j=0; j<sc->params.size_g; j++) {
            index = (leading_pivot + j) % num_p;
            selected[index] += 1;
            sc->gene[i]->pktid[j] = index;
        }
    }
    int coverage = 0;
    for (i=0; i<num_p; i++)
        coverage += selected[i];

    free(selected);
    return coverage;
}

/*
 * Allocate an empty GNC coded packet
 *  gid = -1
 *  coes: zeros
 *  syms: zeros
 */
struct snc_packet *snc_alloc_empty_packet(struct snc_parameters *sp)
{
    struct snc_packet *pkt = calloc(1, sizeof(struct snc_packet));
    if (pkt == NULL)
        return NULL;
    pkt->ucid = -1;
    if (sp->bnc) {
        // For binary code, coding coefficients (bits) are condensed
        pkt->coes = calloc(ALIGN(sp->size_g, 8), sizeof(GF_ELEMENT));
    } else {
        pkt->coes = calloc(sp->size_g, sizeof(GF_ELEMENT));
    }
    if (pkt->coes == NULL)
        goto AllocErr;
    pkt->syms = calloc(sp->size_p, sizeof(GF_ELEMENT));
    if (pkt->syms == NULL)
        goto AllocErr;

    return pkt;

AllocErr:
    snc_free_packet(pkt);
    return NULL;
}

/* Generate a GNC coded packet. Memory is allocated in the function. */
struct snc_packet *snc_generate_packet(struct snc_context *sc)
{
    struct snc_packet *pkt = snc_alloc_empty_packet(&sc->params);
    int gid = schedule_generation(sc);
    encode_packet(sc, gid, pkt);
    return pkt;
}

/*
 * Generate a GNC coded packet in a given memory area.
 * It is the caller's responsibity to allocate memory properly.
 */
int snc_generate_packet_im(struct snc_context *sc, struct snc_packet *pkt)
{
    if (pkt == NULL || pkt->coes == NULL || pkt->syms == NULL)
        return -1;
    if (sc->params.bnc) {
        memset(pkt->coes, 0, ALIGN(sc->params.size_g, 8)*sizeof(GF_ELEMENT));
    } else {
        memset(pkt->coes, 0, sc->params.size_g*sizeof(GF_ELEMENT));
    }
    memset(pkt->syms, 0, sc->params.size_p*sizeof(GF_ELEMENT));
    int gid = schedule_generation(sc);
    encode_packet(sc, gid, pkt);
    return (0);
}

void snc_free_packet(struct snc_packet *pkt)
{
    if (pkt == NULL)
        return;
    if (pkt->coes != NULL)
        free(pkt->coes);
    if (pkt->syms != NULL)
        free(pkt->syms);
    free(pkt);
}


static void encode_packet(struct snc_context *sc, int gid, struct snc_packet *pkt)
{
    pkt->gid = gid;
    int pktid;
    if (sc->params.sys == 1 && sc->count < sc->snum) {
        // send an uncoded source packet
        pktid = sc->count;
        memcpy(pkt->syms, sc->pp[pktid], sc->params.size_p*sizeof(GF_ELEMENT));
        pkt->gid = -1;    // gid=-1 && ucid != -1 indicates it's a systematic packet
        pkt->ucid = pktid;
        sc->nccount[gid] += 1;
        sc->count += 1;
        return;
    }
    /*
    if (sc->params.sys == 1 && sc->nccount[gid] < sc->params.size_b) {
        // Send an uncoded packet
        if (sc->params.bnc == 1) {
            set_bit_in_array(pkt->coes, sc->nccount[gid]);
        } else {
            pkt->coes[sc->nccount[gid]] = 1;
        }
        pktid = sc->gene[gid]->pktid[sc->nccount[gid]];
        memcpy(pkt->syms, sc->pp[pktid], sc->params.size_p*sizeof(GF_ELEMENT));
        pkt->ucid = sc->nccount[gid];      // Mark the uncoded pkt, and store its index amongst the generation
        sc->nccount[gid] += 1;
        return;
    }
    */
    int i;
    GF_ELEMENT co;
    for (i=0; i<sc->params.size_g; i++) {
        pktid = sc->gene[gid]->pktid[i];  // The i-th packet of the gid-th generation
        if (sc->params.bnc) {
            co = (GF_ELEMENT) rand() % 2;                   // Binary network code
            if (co == 1)
                set_bit_in_array(pkt->coes, i);             // Set the corresponding coefficient as 1
        } else {
            co = (GF_ELEMENT) rand() % (1 << 8);     // Randomly generated coding coefficient
            pkt->coes[i] = co;
        }
        galois_multiply_add_region(pkt->syms, sc->pp[pktid], co, sc->params.size_p);
    }
    pkt->ucid = -1;
    sc->nccount[gid] += 1;
    sc->count += 1;
    return;
}

static int schedule_generation(struct snc_context *sc)
{
    if (sc->gnum == 1)
        return 0;

    char *ur = getenv("SNC_NONUNIFORM_RAND");
    if ( ur != NULL && atoi(ur) == 1)
        return banded_nonuniform_sched(sc);
    int gid = rand() % (sc->gnum);
    return gid;
}

/*
 * Non-uniform random scheduling for banded codes
 * NOTE: scheduling of the 0-th and the (M-G)-th generation are not uniform
 * 0-th and (M-G)-th: (G+1)/2M
 * 1-th to (M-G-1)-th: 1/M
 * [G+1, 2, 2, 2,..., 2, G+1]
 * [-----{  2*(M-G-1)  }----]
 */
static int banded_nonuniform_sched(struct snc_context *sc)
{
	int M = sc->snum + sc->cnum;
	int G = sc->params.size_g;
	int upperb = 2*(G+1)+2*(M-G-1);
    int selected = rand() % upperb + 1;
	// int selected = gsl_rng_uniform_int(r, upperb) + 1;

	if (selected <= G+1) {
		selected = 0;
	} else if (selected > (G+1+2*(M-G-1))) {
		selected = sc->gnum - 1;
    } else {
		int residual = selected - (G+1);
		int mod = residual / 2;
		selected = mod + 1;
	}
	return selected;
}

/*
 * Print code summary
 * If called by decoders, it prints overhead and operations as well.
 */
void print_code_summary(struct snc_context *sc, double overhead, double operations)
{
    char typestr[20];
    switch(sc->params.type) {
        case RAND_SNC:
            strcpy(typestr, "RAND");
            break;
        case BAND_SNC:
            strcpy(typestr, "BAND");
            break;
        case WINDWRAP_SNC:
            strcpy(typestr, "WINDWRAP");
            break;
        default:
            strcpy(typestr, "UNKNOWN");
    }
    // precode type
    char *hdpc = getenv("SNC_PRECODE");
    int HDPC = 0;
    if (hdpc != NULL && strcmp(hdpc, "HDPC") == 0)
        HDPC = 1;
    char typestr2[20];
    if (sc->params.size_c == 0) {
        strcpy(typestr2, "NoPrecode");
    } else {
        if (sc->params.bpc) {
            if (HDPC)
                strcpy(typestr2, "BinaryHDPC");
            else
                strcpy(typestr2, "BinaryLDPC");
        } else {
            if (HDPC)
                strcpy(typestr2, "NonBinaryHDPC");
            else
                strcpy(typestr2, "NonBinaryLDPC");
        }
    }
    char typestr3[20];
    if (sc->params.bnc) {
        strcpy(typestr3, "BinaryNC");
    } else {
        strcpy(typestr3, "NonBinaryNC");
    }
    char typestr4[20];
    if (sc->params.sys) {
        strcpy(typestr4, "Systematic");
    } else {
        strcpy(typestr4, "NonSystematic");
    }
    printf("datasize: %d ", sc->params.datasize);
    printf("size_p: %d ", sc->params.size_p);
    printf("snum: %d ", sc->snum);
    printf("size_c: %d ", sc->params.size_c);
    printf("size_b: %d ", sc->params.size_b);
    printf("size_g: %d ", sc->params.size_g);
    printf("type: [%s::%s::%s::%s] ", typestr, typestr2, typestr3, typestr4);
    printf("gnum: %d ", sc->gnum);
    if (operations != 0) {
        printf("overhead: %.6f ", overhead);
        printf("computation: %.4f\n", operations);
    } else {
        printf("\n");
    }
}

