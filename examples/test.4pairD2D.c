#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "Simulate RLNC network coded D2D cooperative transmissions\n\
                \n\
                \n\
                    --------U1\n\
                   /        ||\n\
                S -         ||\n\
                   \\       ||\n\
                    --------U2\n\
                \n\
                usage: GF_SIZE=q ./programName pktsize snum alpha\n";

int main(int argc, char *argv[])
{
    if (argc != 4) {
        printf("%s\n", usage);
        exit(1);
    }
    int sched_t = MLPI_SCHED;

    struct snc_parameters sp;
    sp.type = BAND_SNC;
    int decoder_t = CBD_DECODER;
    sp.datasize = atoi(argv[1]) * atoi(argv[2]);
    sp.size_p   = atoi(argv[1]);
    sp.size_c   = 0;
    sp.size_b   = atoi(argv[2]);
    sp.size_g   = atoi(argv[2]);
    sp.bpc      = 0;
    char *gf_evar = getenv("GF_SIZE");
    if ( gf_evar != NULL && atoi(gf_evar) == 2) {
        sp.bnc = 1;
    } else {
        sp.bnc = 0;
    }
    sp.sys      = 1;
    sp.seed     = -1;

    double alpha = atof(argv[3]);
	double delta[8] = {0.02, 0.2, 0.3, 0.01, 0.05, 0.2, 0.05, 0.2};  // broadcast link loss probability
	double epsil[8] = {0.01, 0.01, 0.05, 0.01, 0.02, 0.01, 0.01, 0.02};  // broadcast link loss probability
    double pD2D[4] = {0.2, 0.3, 0.2, 0.2};
    int nD2D[8] = {1, 0, 0, 1, 1, 0, 1, 0};
    int nD2Dmax[4];
    nD2Dmax[0] = nD2D[0] > nD2D[1] ? nD2D[0] : nD2D[1];
    nD2Dmax[1] = nD2D[2] > nD2D[3] ? nD2D[2] : nD2D[3];
    nD2Dmax[2] = nD2D[4] > nD2D[5] ? nD2D[4] : nD2D[5];
    nD2Dmax[3] = nD2D[6] > nD2D[7] ? nD2D[6] : nD2D[7];
	
    int i, j;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    char *buf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    /* Create GNC encoding context */
    if ((sc = snc_create_enc_context(buf, &sp)) == NULL) {
        fprintf(stderr, "Cannot create snc_context.\n");
        exit(1);
    }

    int nusers = 2;
    int npairs = 4;
    // Create decoders and buffer for each user of each pair
    struct snc_decoder **decoders =  malloc(sizeof(struct snc_decoder*) * nusers * npairs);
    struct snc_buffer **buffers = malloc(sizeof(struct snc_buffer*) * nusers * npairs);
    for (i=0; i<nusers*npairs; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
        buffers[i] = snc_create_buffer(&sp, sp.size_g);  // create a large enough buffer
    }

    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nusers*npairs);  // count network uses of users
    double *Ecoop = calloc(sizeof(double), nusers*npairs);  // record cooperation energy consumption of each user

    int decoded = 0;
    int count = 0;
    while (decoded != nusers*npairs) {
        count++;
        pkt = snc_generate_packet(sc);

        // Phase I: broadcast to n users
        for (i=0; i<nusers*npairs; i++) {
            if (snc_decoder_finished(decoders[i])) 
                continue;  // skip completed users

            nuse[i] += 1;
            if (rand() % 100 >= delta[i] * 100) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                struct snc_packet *buf_pkt = snc_duplicate_packet(pkt, &sp);
                snc_buffer_packet(buffers[i], buf_pkt);              // save a copy in the buffer, for the purpose of D2D comm.
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
        
        // Phase II: D2D communication
        // D2D phase occurs with probability pD2D, and if occurred, each user sends a packet to the other
        for (int k=0; k<npairs; k++) {
            if (pD2D[k] != 0 && rand() % 100 < pD2D[k]*100) {
                // Cooperation among each pair
                for (int n=0; n<nD2Dmax[n]; n++) {
                    // U1 send to U2
                    if (n<nD2D[k*nusers] && !snc_decoder_finished(decoders[k*nusers+1])) {
                        pkt = snc_recode_packet(buffers[k*nusers], sched_t);
                        if (pkt == NULL)
                            continue;
                        Ecoop[k*nusers] += alpha;    // counting energy consumption for the cooperation packet
                        if (rand() % 100 >= epsil[k*nusers]) {
                            struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                            snc_process_packet(decoders[k*nusers+1], d2d_pkt);
                            snc_free_packet(d2d_pkt);
                        }
                        snc_free_packet(pkt);
                    }

                    // U2 send to U1
                    if (n<nD2D[k*nusers+1] && !snc_decoder_finished(decoders[k*nusers])) {
                        pkt = snc_recode_packet(buffers[k*nusers+1], sched_t);
                        if (pkt == NULL)
                            continue;
                        Ecoop[k*nusers+1] += alpha;    // counting energy consumption for the cooperation packet
                        if (rand() % 100 >= epsil[k*nusers+1]) {
                            struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                            snc_process_packet(decoders[k*nusers], d2d_pkt);
                            snc_free_packet(d2d_pkt);
                        }
                        snc_free_packet(pkt);
                    }
                }
            }
        }

        // check decoder status
        decoded = 0;
        for (i=0; i<npairs*nusers; i++) {
            if (snc_decoder_finished(decoders[i]))
                decoded += 1;
        }
    }

    for (i=0; i<npairs*nusers; i++) {
        struct snc_context *dsc = snc_get_enc_context(decoders[i]);
        unsigned char *rec_buf = snc_recover_data(dsc);
        if (memcmp(buf, rec_buf, sp.datasize) != 0)
            fprintf(stderr, "recovered is NOT identical to original.\n");
        printf("user: %d overhead: %.6f ops: %.6f network-uses: %d\n", i, snc_decode_overhead(decoders[i]), snc_decode_cost(decoders[i]), nuse[i]);
    }
    int GF_size = snc_get_GF_size(&sp);
    printf("snum: %d pktsize: %d GF_SIZE: %d nBroadcast: %d CoopEnergy: %.6f \n", sp.datasize/sp.size_p, sp.size_p, GF_size, count);
    for (i=0; i<npairs; i++) {
        printf("Pair-%d: ", i);
        for (j=0; j<nusers; j++)
            printf(" delta[%d]: %.3f ", j+1, delta[i*nusers+j]);
        for (j=0; j<nusers; j++)
            printf(" epsil[%d]: %.3f ", j+1, epsil[i*nusers+j]);
        printf(" pD2D: %.3f ", pD2D[i]);
        for (j=0; j<nusers; j++)
            printf(" nD2D[%d]: %d ", j+1, nD2D[i*nusers+j]);
        for (j=0; j<nusers; j++)
            printf(" nuse[%d]: %d ", j+1, nuse[i*nusers+j]);
        printf(" CoopEnergy: %.6f ", Ecoop[i*nusers]+Ecoop[i*nusers+1]);
        printf("\n");
    }
    snc_free_enc_context(sc);
    for (i=0; i<nusers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}

