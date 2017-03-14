#include <math.h>
#include "common.h"
#include "galois.h"
#include "sparsenc.h"

/* Schedule a subgeneration to recode a packet according
 * to the specified scheduling type. */
static int schedule_recode_generation(struct snc_buffer *buf, int sched_t);
static int banded_nonuniform_sched(struct snc_buffer *buf);

// Needed for buffering systematic codes
static struct snc_context *sc;          // encoding context duplicated at the recoder
static ID_list **gene_nbr = NULL;       // Lists of subgen neighbors of each packet
                                        // Needed if snc is systematic

struct snc_buffer *snc_create_buffer(struct snc_parameters *sp, int bufsize)
{
    static char fname[] = "snc_create_decoder";
    int i;
    struct snc_buffer *buf;
    if ((buf = calloc(1, sizeof(struct snc_buffer))) == NULL) {
        fprintf(stderr, "calloc snc_buffer fail\n");
        goto Error;
    }
    //initialize recoding buffer
    buf->params = *sp;
    // determine number of generations with sp
    int num_src = ALIGN(buf->params.datasize, buf->params.size_p);
    buf->snum = num_src;
    int num_chk = sp->size_c;
    buf->cnum = sp->size_c;
    if (buf->params.type == BAND_SNC)
        buf->gnum  = ALIGN((num_src+num_chk-buf->params.size_g), buf->params.size_b) + 1;
    else
        buf->gnum  = ALIGN( (num_src+num_chk), buf->params.size_b);

    buf->size = bufsize;
    buf->nemp = 0;
    if ((buf->gbuf = calloc(buf->gnum, sizeof(struct snc_packet **))) == NULL) {
        fprintf(stderr, "%s: calloc buf->gbuf\n", fname);
        goto Error;
    }
    for (i=0; i<buf->gnum; i++) {
        /* Initialize pointers of buffered packets of each generation as NULL */
        if ((buf->gbuf[i] = calloc(bufsize, sizeof(struct snc_packet *))) == NULL) {
            fprintf(stderr, "%s: calloc buf->gbuf[%d]\n", fname, i);
            goto Error;
        }
    }
    if ((buf->nc = calloc(buf->gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->nc\n", fname);
        goto Error;
    }
    if ((buf->pn = calloc(buf->gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->pn\n", fname);
        goto Error;
    }
    if ((buf->nsched = calloc(buf->gnum, sizeof(int))) == NULL) {
        fprintf(stderr, "%s: calloc buf->nsched\n", fname);
        goto Error;
    }
    if (sp->sys == 1) {
        // allocate buffer for uncoded packets
        if ((buf->sysbuf = calloc(buf->snum, sizeof(struct snc_packet *))) == NULL){
            fprintf(stderr, "%s: calloc buf->sysbuf\n", fname);
            goto Error;
        }
        buf->sysnum = 0;
        buf->sysptr = 0;
        sc = snc_create_enc_context(NULL, sp);  //encoding context is needed for systematic forwarding/recoding
        gene_nbr = build_subgen_nbr_list(sc);
    }
    return buf;

Error:
    snc_free_buffer(buf);
    return NULL;
}

/*
 * Buffer structure example (nc=1, pn=1)
 * snc_packet
 *      ^           NULL         NULL
 *      |            |            |
 *      |            |            |
 * gbuf[gid][0] gbuf[gid][1] gbuf[gid][2] ......... gbuf[gid][size-1]
 *                   ^
 *                   |
 *                   |
 *                pn = 1
 */
void snc_buffer_packet(struct snc_buffer *buf, struct snc_packet *pkt)
{
    int gid = pkt->gid;
    if (gid == -1 && pkt->ucid != -1) {
        // This is a systematic packet
        buf->sysbuf[buf->sysnum] = pkt;
        buf->sysnum += 1;
        // Count the packet as a received packet for each subgen containing it
        ID *item = gene_nbr[pkt->ucid]->first;
        while (item != NULL) {
            gid = item->data;
            if (buf->nc[gid] == 0)
                buf->nemp++;
            buf->nc[gid]++;
            item = item->next;
        }
        return;
    }
    if (buf->nc[gid] == 0) {
        // Buffer of the generation is empty
        buf->gbuf[gid][0] = pkt;
        buf->nc[gid]++;
        buf->nemp++;
    } else if (buf->nc[gid] == buf->size) {
        // Buffer of the generation is full, remove in a FIFO manner
        snc_free_packet(buf->gbuf[gid][buf->pn[gid]]);  //discard packet previously stored in the position
        buf->gbuf[gid][buf->pn[gid]] = pkt;
    } else {
        // Buffer is neither empty nor full
        buf->gbuf[gid][buf->pn[gid]] = pkt;
        buf->nc[gid]++;
    }
    // Update position for next incoming coded packets
    buf->pn[gid] = (buf->pn[gid] + 1) % buf->size;
    return;
}

struct snc_packet *snc_recode_packet(struct snc_buffer *buf, int sched_t)
{
    struct snc_packet *pkt = snc_alloc_empty_packet(&buf->params);
    if (pkt == NULL)
        return NULL;

    if (snc_recode_packet_im(buf, pkt, sched_t) == 0)
        return pkt;
    else
        return NULL;
    return pkt;
}

int snc_recode_packet_im(struct snc_buffer *buf, struct snc_packet *pkt, int sched_t)
{

    int gid = schedule_recode_generation(buf, sched_t);
    if (gid == -1)
        return -1;

    if (gid == buf->gnum) {
        // A systematic packet is scheduled
        pkt->gid = -1;
        pkt->ucid = buf->sysbuf[buf->sys_sched]->ucid;
        memset(pkt->syms, 0, sizeof(GF_ELEMENT)*buf->params.size_p);
        memcpy(pkt->syms, buf->sysbuf[buf->sys_sched]->syms, sizeof(GF_ELEMENT)*buf->params.size_p);
        return 0;
    }

    // Generate a normal recoded GNC packet
    pkt->gid = gid;
    pkt->ucid = -1;
    // Clean up pkt
    if (buf->params.bnc) {
        memset(pkt->coes, 0, ALIGN(buf->params.size_g, 8)*sizeof(GF_ELEMENT));
    } else {
        memset(pkt->coes, 0, buf->params.size_g*sizeof(GF_ELEMENT));
    }
    memset(pkt->syms, 0, sizeof(GF_ELEMENT)*buf->params.size_p);
    GF_ELEMENT co = 0;
    int i;
    // First, go through systematic packet list, combine those belonging to the shceduled generation
    int ucnum = 0;  // num of systematic packets belonging to this gid
    for (i=0; i<buf->sysnum; i++) {
        if (!exist_in_list(gene_nbr[buf->sysbuf[i]->ucid], gid))
            continue;            // skip packets not belong to gid
        // Find the sys packet's corresponding index in the generation.
        // If buf->sysbuf[i]->ucid doesn't belong to the generation, skip.
        int relative_idx = has_item(sc->gene[gid]->pktid, buf->sysbuf[i]->ucid, sc->params.size_g); 
        if (relative_idx == -1)
            continue;
        ucnum += 1;
        if (sc->params.bnc) {
            co = (GF_ELEMENT) rand() % 2;                   // Binary network code
            if (co == 1)
                set_bit_in_array(pkt->coes, relative_idx);  // Set the corresponding coefficient as 1
        } else {
            co = (GF_ELEMENT) rand() % (1 << 8);     // Randomly generated coding coefficient
            pkt->coes[relative_idx] = co;
        }
        galois_multiply_add_region(pkt->syms, buf->sysbuf[i]->syms, co, sc->params.size_p);
    }

    // Second, go through the buffered coded packets of the generation
    // Note that we need to subtract ucnum from buf->nc[gid]
    for (i=0; i<buf->nc[gid]-ucnum; i++) {
        if (buf->params.bnc == 1) {
            co = rand() % 2;
            if (co == 1)
                galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, ALIGN(buf->params.size_g, 8));
        } else {
            co = rand() % (1 << 8);
            galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, buf->params.size_g);
        }
        galois_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, buf->params.size_p);
    }
    return 0;
}

void snc_free_buffer(struct snc_buffer *buf)
{
    if (buf == NULL)
        return;
    int i;
    for (i=0; i<buf->gnum; i++) {
        if (buf->gbuf != NULL && buf->gbuf[i] != NULL) {
            /* Free bufferred packets, if any */
            for (int j=0; j<buf->size; j++)
                snc_free_packet(buf->gbuf[i][j]);
            /* Free the pointer array */
            free(buf->gbuf[i]);
        }
    }
    if (buf->gbuf != NULL)
        free(buf->gbuf);
    if (buf->nc != NULL)
        free(buf->nc);
    if (buf->pn != NULL)
        free(buf->pn);
    if (buf->nsched != NULL)
        free(buf->nsched);
    // Free systematic packets
    if (buf->sysbuf != NULL) {
        for (i=0; i<buf->snum; i++) {
            // Free systematic packets
            if (buf->sysbuf[i] != NULL)
                snc_free_packet(buf->sysbuf[i]);
        }
        free(buf->sysbuf);
    }
    free(buf);
    buf = NULL;

    // free static variables if used
    free_subgen_nbr_list(sc, gene_nbr);
    snc_free_enc_context(sc);
    return;
}

/*
 *
 * Return:
 *  -1            - scheduling failed
 *  [0, gnum-1]   - scheduled generation
 *  gnum          - forward systematic packet
 */
static int schedule_recode_generation(struct snc_buffer *buf, int sched_t)
{

    if ((sched_t == RAND_SCHED_SYS || sched_t == MLPI_SCHED_SYS) && (buf->params.sys != 1)) {
        fprintf(stderr, "Systematic scheduling can only be used with systematic coding schemes.\n");
        return -1;
    }

    if (buf->nemp == 0 && buf->sysnum == 0)
        return -1;

    int gid;

    if (sched_t == RAND_SCHED_SYS || sched_t == MLPI_SCHED_SYS) {
        if (buf->sysptr < buf->sysnum) {
            buf->sys_sched = buf->sysptr;
            buf->sysptr += 1;
            // update nsched of subgens containing the packet
            ID *item = gene_nbr[buf->sysbuf[buf->sys_sched]->ucid]->first;
            while (item != NULL) {
                gid = item->data;
                buf->nsched[gid] += 1;
                item = item->next;
            }
            return buf->gnum;
        }
    }

    if (sched_t == TRIV_SCHED) {
        gid = rand() % buf->gnum;
        buf->nsched[gid]++;
        return gid;
    }

    if (sched_t == RAND_SCHED || sched_t == RAND_SCHED_SYS) {
        int index = rand() % buf->nemp;
        int i = -1;
        gid = 0;
        while ( i != index) {
            if (buf->nc[gid++] != 0)
                i++;
        }
        buf->nsched[gid-1]++;
        return gid-1;
    }

    if (sched_t == MLPI_SCHED || sched_t == MLPI_SCHED_SYS) {
        gid = 0;
        int max = buf->nc[gid] - buf->nsched[gid];
        for (int j=0; j<buf->gnum; j++) {
            if (buf->nc[j] - buf->nsched[j] > max) {
                max = buf->nc[j] - buf->nsched[j];
                gid = j;
            }
        }
        buf->nsched[gid]++;
        return gid;
    }

    if (sched_t == NURAND_SCHED) {
        return banded_nonuniform_sched(buf);
    }
}


/*
 * Non-uniform random scheduling for banded codes
 * NOTE: scheduling of the 0-th and the (M-G)-th generation are not uniform
 * 0-th and (M-G)-th: (G+1)/2M
 * 1-th to (M-G-1)-th: 1/M
 * [G+1, 2, 2, 2,..., 2, G+1]
 * [-----{  2*(M-G-1)  }----]
 */
static int banded_nonuniform_sched(struct snc_buffer *buf)
{
	int M = buf->snum + buf->cnum;
	int G = buf->params.size_g;
	int upperb = 2*(G+1)+2*(M-G-1);

	int found = 0;
	int selected = -1;
	while (found ==0) {
		selected = (rand() % upperb) + 1;

		if (selected <= G+1) {
			selected = 0;
		} else if (selected > (G+1+2*(M-G-1))){
			selected = buf->gnum - 1;
		} else {
			int residual = selected - (G+1);
			int mod = residual / 2;
			selected = mod + 1;
		}
		if (buf->nc[selected] != 0)
			found = 1;
	}
	return selected;
}
