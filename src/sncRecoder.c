#include <math.h>
#include "common.h"
#include "galois.h"
#include "sparsenc.h"

/* Schedule a subgeneration to recode a packet according
 * to the specified scheduling type. */
static int schedule_recode_generation(struct snc_buffer *buf, int sched_t);
static int banded_nonuniform_sched(struct snc_buffer *buf);

// Needed for buffering systematic codes
static struct snc_context *sc = NULL;          // encoding context duplicated at the recoder
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
        buf->newsys = -1;
        buf->sysgid = -1;
        buf->sysidx = -1;
        /*
        // allocate buffer for uncoded packets
        if ((buf->sysbuf = calloc(buf->size, sizeof(struct snc_packet *))) == NULL){
            fprintf(stderr, "%s: calloc buf->sysbuf\n", fname);
            goto Error;
        }
        buf->spn    = 0;
        buf->sysnum = 0;
        buf->newsys = -1;
        buf->sysptr = 0;
        */
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
    int gfpower = buf->params.gfpower;
    int i, j;
    int gid = pkt->gid;
    // Receive and buffer a systematic packet
    if (gid == -1 && pkt->ucid != -1) {
        // find which subgens it belongs to 
        ID *item = gene_nbr[pkt->ucid]->first;
        int sgid;
        int set_sched = 0;
        while (item != NULL) {
            sgid = item->data;
            // Duplicate and store in the buffers of the subgen to which the packet belongs
            // i.e., treat is as a normal coded packet
            struct snc_packet *pktcopy = snc_duplicate_packet(pkt, &buf->params);
            int relative_idx = has_item(sc->gene[sgid]->pktid, pkt->ucid, sc->params.size_g);
            // We need to set the corresponding coding coefficient to 1
            memset(pktcopy->coes, 0, ALIGN(buf->params.size_g * gfpower, 8) * sizeof(GF_ELEMENT));
            if (buf->params.gfpower == 1) {
                set_bit_in_array(pktcopy->coes, relative_idx);
            } else if (buf->params.gfpower == 8){
                pktcopy->coes[relative_idx] = 1;
            } else {
                pack_bits_in_byte_array(pktcopy->coes, ALIGN(buf->params.size_g*gfpower, 8), 1, gfpower, relative_idx);
            }
            if (buf->gbuf[sgid][buf->pn[sgid]] != NULL) {
                snc_free_packet(buf->gbuf[sgid][buf->pn[sgid]]);  //discard packet previously stored in the position
            }
            buf->gbuf[sgid][buf->pn[sgid]] = pktcopy;
            if (set_sched == 0) {
                // Update information in case of systematic scheduling
                // Here we pick the first subgen (it doesn't matter which one to choose)
                buf->newsys = 1;
                buf->sysgid = sgid;
                buf->sysidx = buf->pn[sgid];
                set_sched = 1;
            }
            if (buf->nc[sgid] == 0) {
                buf->nemp++;
            }
            if (buf->nc[sgid]<buf->size) {
                buf->nc[sgid]++;
            }

            // Update position for next incoming packet of the subgen
            buf->pn[sgid] = (buf->pn[sgid] + 1) % buf->size;
            
            item = item->next;
        }
        snc_free_packet(pkt);
        return;
    }

    // Receive and buffer a normal coded packet
    // Reset newsys immediately: we would only schedule a systematic packet if it is the latest
    // received packet.
    buf->newsys = -1;
    buf->sysgid = -1;
    buf->sysidx = -1;
    if (buf->nc[gid] == 0) {
        // Buffer of the generation is empty
        buf->gbuf[gid][0] = pkt;
        buf->nc[gid]++;
        buf->nemp++;
    } else if (buf->nc[gid] < buf->size) {
        // the buffer is not empty nor full
        buf->gbuf[gid][buf->pn[gid]] = pkt;
        buf->nc[gid]++;
    } else {
        // The buffer is full
        // Strategy (a): buffer in a FIFO manner
        /*
        snc_free_packet(buf->gbuf[gid][buf->pn[gid]]);  //discard packet previously stored in the position
        buf->gbuf[gid][buf->pn[gid]] = pkt;
        */
        // Strategy (b): buffer in the accumulator manner, i.e., randomly add the packet to the
        // buffered packets. Ref. Lun2006 "An analysis of finite-memory random linear coding on packet streams"
        // printf("Buffer full, accumulating the received packet to buffered packets\n");
         
        for (int i=0; i<buf->size; i++) {
            GF_ELEMENT co = (GF_ELEMENT) genrand_int32() % (1<<gfpower);
            if (gfpower == 1) {
                if (co == 1) {
                    galois_multiply_add_region(buf->gbuf[gid][i]->coes, pkt->coes, co, ALIGN(buf->params.size_g, 8));
                }
            } else if (gfpower == 8) {
                galois_multiply_add_region(buf->gbuf[gid][i]->coes, pkt->coes, co, buf->params.size_g);
            } else {
                // coding coefficients of 2,3,...,7 bits
                galois2n_multiply_add_region(buf->gbuf[gid][i]->coes, pkt->coes, co, gfpower, buf->params.size_g, ALIGN(buf->params.size_g*gfpower, 8));
            }

            // multiply_add the coded symbols
            if (gfpower == 1 || gfpower == 8) {
                galois_multiply_add_region(buf->gbuf[gid][i]->syms, pkt->syms, co, buf->params.size_p);
            } else {
                int nelem = ALIGN(buf->params.size_p*8, gfpower);
                galois2n_multiply_add_region(buf->gbuf[gid][i]->syms, pkt->syms, co, gfpower, nelem, buf->params.size_p);
            }
        }
        
        
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
    else {
        snc_free_packet(pkt);
        return NULL;
    }
}

int snc_recode_packet_im(struct snc_buffer *buf, struct snc_packet *pkt, int sched_t)
{
    int gfpower = buf->params.gfpower;
    int gid = schedule_recode_generation(buf, sched_t);
    if (gid == -1)
        return -1;
    
    // Clean up pkt
    memset(pkt->coes, 0, ALIGN(buf->params.size_g * gfpower, 8)*sizeof(GF_ELEMENT));
    memset(pkt->syms, 0, sizeof(GF_ELEMENT)*buf->params.size_p);
    
    if (gid == buf->gnum) {
        // A systematic packet is scheduled
        pkt->gid = -1;
        pkt->ucid = buf->gbuf[buf->sysgid][buf->sysidx]->ucid;
        // memcpy(pkt->coes, buf->gbuf[buf->sysgid][buf->sysidx]->coes, ALIGN(buf->params.size_g*gfpower, 8)*sizeof(GF_ELEMENT));
        memcpy(pkt->syms, buf->gbuf[buf->sysgid][buf->sysidx]->syms, buf->params.size_p*sizeof(GF_ELEMENT));
        buf->newsys = -1;
        buf->sysgid = -1;
        buf->sysidx = -1;
        return 0;
    }

    // Generate a normal recoded GNC packet
    pkt->gid = gid;
    pkt->ucid = -1;
    GF_ELEMENT co = 0;
    int i, j;
    // Go through the buffered packets of the subgeneration
    for (i=0; i<buf->nc[gid]; i++) {
        co = (GF_ELEMENT) genrand_int32() % (1 << gfpower);
        if (co == 0) {
            continue;
        }
        if (gfpower == 1) {
            galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, ALIGN(buf->params.size_g, 8));
        } else if (gfpower == 8) {
            galois_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, buf->params.size_g);
        } else {
            // coding coefficients of 2,3,...,7 bits
            // int nelem = ALIGN(buf->params.size_g*8, gfpower);
            galois2n_multiply_add_region(pkt->coes, buf->gbuf[gid][i]->coes, co, gfpower, buf->params.size_g, ALIGN(buf->params.size_g*gfpower, 8));
        }
        // linear combinations of coded symbols
        if (gfpower == 1 || gfpower == 8) {
            galois_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, buf->params.size_p);
        } else {
            int nelem = ALIGN(buf->params.size_p*8, gfpower);
            galois2n_multiply_add_region(pkt->syms, buf->gbuf[gid][i]->syms, co, gfpower, nelem, buf->params.size_p);
        }
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
    free(buf);
    buf = NULL;

    // free static variables if used
    if (gene_nbr != NULL && sc != NULL) {
        free_subgen_nbr_list(sc, gene_nbr);
        gene_nbr = NULL;
        snc_free_enc_context(sc);
        sc = NULL;
    }
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

    /*
     * If code is not systematic but sched_t was to systematic, it will automatically
     * fall back to non-systematic scheduling.
    if ((sched_t == RAND_SCHED_SYS || sched_t == MLPI_SCHED_SYS) && (buf->params.sys != 1)) {
        fprintf(stderr, "Systematic scheduling can only be used with systematic coding schemes.\n");
        return -1;
    }
    */

    if (buf->nemp == 0)
        return -1;

    int gid;

    if (sched_t == RAND_SCHED_SYS || sched_t == MLPI_SCHED_SYS) {
        if (buf->newsys != -1) {
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
