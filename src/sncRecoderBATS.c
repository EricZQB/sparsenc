#include "common.h"
#include "galois.h"
#include "sparsenc.h"

static int s_neq_r = 0;  // indicate whether the sending and receiving batch do not match.
                     // i.e., whether there are more than 1 batches buffered
static int s_count = 0;  // count num of pkts sent from the current sending batch

void visualize_buffer2(struct snc_buffer_bats *buf);

struct snc_buffer_bats *snc_create_buffer_bats(struct snc_parameters *sp, int bufsize)
{
    static char fname[] = "snc_create_buffer_bats";

    struct snc_buffer_bats *buf = malloc(sizeof(struct snc_buffer_bats));

    buf->params = *sp;
    buf->srbuf = calloc(bufsize, sizeof(struct snc_packet*));  // allocate buffer packet pointers
    buf->bufsize = bufsize;
    buf->sbatchid = -1;  // empty buffer
    buf->s_first = -1;
    buf->r_last = -1;
    return buf;
}

void snc_buffer_packet_bats(struct snc_buffer_bats *buf, struct snc_packet *pkt)
{
    int pos = -1;   // Pos index where new packet is stored
    //printf("buffering batch: %d\n", pkt->batchid);
    // the very first batch received
    if (buf->sbatchid == -1) { 
        pos = 0;
        buf->r_last  = pos;
        buf->s_first = pos;
        buf->sbatchid = pkt->gid;
        buf->srbuf[pos] = pkt;
        return;
    }

    if (pkt->gid != buf->sbatchid)
        s_neq_r = 1;

    if (((buf->r_last+1) % buf->bufsize) == buf->s_first) {
        // buffer is full, replace the oldest packet with the new arrived one
        int new_sfirst = (buf->s_first+1) % buf->bufsize;
        if (buf->srbuf[new_sfirst]->gid != buf->sbatchid) {
            printf("\n");
            visualize_buffer2(buf);
        }
        pos = buf->s_first;
        snc_free_packet(buf->srbuf[pos]);
        buf->srbuf[pos] = NULL;
        buf->srbuf[pos] = pkt;
        buf->r_last = pos;
        
        if (buf->srbuf[new_sfirst]->gid != buf->sbatchid) {
            printf("%d packets sent from batch %d ", s_count, buf->sbatchid);
            buf->s_first = new_sfirst;
            buf->sbatchid = buf->srbuf[buf->s_first]->gid; // sending batch id is updated to the new s_first packet
            printf("change send batch to %d [full-buffer]\n", buf->sbatchid);
            s_count = 0;
            s_neq_r = 0;
            visualize_buffer2(buf);
            printf("\n");
        } else {
            buf->s_first = new_sfirst;   // update s_first anyway
        }
        return;
    } else {
        pos = (buf->r_last + 1) % buf->bufsize;
        // just double check
        if (buf->srbuf[pos] != NULL) {
            snc_free_packet(buf->srbuf[pos]);
            buf->srbuf[pos] = NULL;
        }
        buf->srbuf[pos] = pkt;
        buf->r_last = pos;
        return;
    }
    

    return ;
}

struct snc_packet *snc_recode_packet_bats(struct snc_buffer_bats *buf)
{
    struct snc_packet *pkt = snc_alloc_empty_packet(&buf->params);
    if (pkt == NULL)
        return NULL;

    if (snc_recode_packet_bats_im(buf, pkt) == 0)
        return pkt;
    else {
        snc_free_packet(pkt);
        return NULL;
    }
}

// Recode a BATS coded packet in a pre-allocated snc_packet space
int snc_recode_packet_bats_im(struct snc_buffer_bats *buf, struct snc_packet *pkt)
{
    if (buf->sbatchid < 0) {
        //printf("Buffer has no batch buffered yet\n");
        return -1;
    }

    int s_pos = buf->s_first; // first packet of the sending batch
    int i, pos;
    
    if (s_count >= buf->params.size_b && s_neq_r == 0) {
        printf("I have sent %d packets for batch %d, not going to send more.\n", s_count, buf->sbatchid);
        // There is only one batch buffered, but it has sent BTS packets.
        return -1;      // This is to follow the recoding algorithm strictly. Do not recode
                          // is sub-optimal since we're wasting transmission opportunity. But 
                          // I still do this just want to see the lower bound of the performance.
    }
    // change sending batch if the curren sending batch has sent BTS packets
    if (s_count >= buf->params.size_b && s_neq_r == 1) {
        printf("\n");
        visualize_buffer2(buf);
        for (i=0; i<buf->bufsize; i++) {
            pos = (s_pos + i) % buf->bufsize;
            if (buf->srbuf[pos] != NULL && buf->srbuf[pos]->gid != buf->sbatchid) {
                printf("%d packets sent from batch %d, change sending batch to %d [max-sent]\n", s_count, buf->sbatchid, buf->srbuf[pos]->gid);
                buf->sbatchid = buf->srbuf[pos]->gid;
                buf->s_first = pos;
                s_pos = pos;
                s_count = 0;
                break;
            }
        }
        visualize_buffer2(buf);
        printf("\n");
        s_neq_r = 0; // reset s_neqr_r. If there are >2 batches in the buffer,
                     // it will be set again when next packet is received. 
    }

    // Reset the snc_packet structure
    pkt->gid = buf->sbatchid;
    pkt->ucid = -1;
    // Clean up pkt
    if (buf->params.bnc) {
        memset(pkt->coes, 0, ALIGN(buf->params.size_g, 8)*sizeof(GF_ELEMENT));
    } else {
        memset(pkt->coes, 0, buf->params.size_g*sizeof(GF_ELEMENT));
    }
    memset(pkt->syms, 0, sizeof(GF_ELEMENT)*buf->params.size_p);
    // Recoding
    GF_ELEMENT co =0;    
    for (i=0; i<buf->bufsize; i++) {
        pos = (s_pos + i) % buf->bufsize;
        if (buf->srbuf[pos] == NULL || buf->srbuf[pos]->gid != buf->sbatchid)
            break;      // packets belonging to the same batch must be stored adjacently (TODO: is this really true in the asynchronous mode?).
        // Perform random linear combination of buffered packets belonging to the same batch
        if (buf->params.bnc == 1) {
            if ((co = rand() % 2) == 1)
                galois_multiply_add_region(pkt->coes, buf->srbuf[pos]->coes, co, ALIGN(buf->params.size_g, 8));
        } else {
            co = rand() % (1 << 8);
            galois_multiply_add_region(pkt->coes, buf->srbuf[pos]->coes, co, buf->params.size_g);
        }
        galois_multiply_add_region(pkt->syms, buf->srbuf[pos]->syms, co, buf->params.size_p);
    }
    s_count += 1;
    return 0;
}


void snc_free_buffer_bats(struct snc_buffer_bats *buf)
{
    if (buf == NULL)
        return;
    int i;
    for (i=0; i<buf->bufsize; i++) {
        if (buf->srbuf[i] != NULL)
            snc_free_packet(buf->srbuf[i]);
    }
    free(buf->srbuf);
    free(buf);
    buf = NULL;
    return;
}

void visualize_buffer2(struct snc_buffer_bats *buf)
{
    printf("buffer size: %d sbatchid: %d s_first: %d r_last: %d s_count: %d s_neq_r: %d\n", buf->bufsize, buf->sbatchid, buf->s_first, buf->r_last, s_count, s_neq_r);
    for (int i=0; i<buf->bufsize; i++) {
        if (buf->srbuf[i] != NULL) {
            printf("%d\t", buf->srbuf[i]->gid);
        } else {
            printf("-1\t");
        }
    }
    printf("\n");
}
