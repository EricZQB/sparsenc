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
                usage: ./programName nuser broadcast_link_loss_probs cooperation_params\n";

int main(int argc, char *argv[])
{
    if (argc < 3) {
        printf("%s\n", usage);
        exit(1);
    }
    int sched_t = MLPI_SCHED;


    struct snc_parameters sp;
    sp.type = BAND_SNC;
    int decoder_t = CBD_DECODER;
    int M = 16;
    int pktsize = 200;
    sp.datasize = M * pktsize;
    sp.size_p   = pktsize;
    sp.size_c   = 0;
    sp.size_b   = M;
    sp.size_g   = M;
    sp.bpc      = 0;
    sp.gfpower  = 8;
    sp.sys      = 1;
    sp.seed     = -1;

    double alpha = 0.05;
    int nUsers = atoi(argv[1]);
    double *delta = calloc(nUsers, sizeof(double)); // broadcast link loss probability
    double *pD2D  = calloc(nUsers/2, sizeof(double));
    int    *nD2D  = calloc(nUsers, sizeof(int));
    double *epsil = calloc(nUsers/2, sizeof(double)); // loss probabilities from the head user to each of the others
    int i, j;
    for (i=0; i<nUsers; i++) {
        delta[i] = atof(argv[2+i]);
    }
    int *nD2Dmax = calloc(nUsers/2, sizeof(int));      // max n of each pair
    for (i=0; i<nUsers/2; i++) {
        pD2D[i] = atof(argv[nUsers+2+i*4]);
        nD2D[i*2] = atoi(argv[nUsers+2+i*4+1]);
        nD2D[i*2+1] = atoi(argv[nUsers+2+i*4+2]);
        epsil[i] = atof(argv[nUsers+2+i*4+3]);
        nD2Dmax[i] = nD2D[i*2] > nD2D[i*2+1] ? nD2D[i*2] : nD2D[i*2+1];
    }
	
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

    int npairs = nUsers/2;
    // Create decoders and buffer for each user of each pair
    struct snc_decoder **decoders =  malloc(sizeof(struct snc_decoder*) * nUsers);
    struct snc_buffer **buffers = malloc(sizeof(struct snc_buffer*) * nUsers);
    for (i=0; i<nUsers; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
        buffers[i] = snc_create_buffer(&sp, sp.size_g);  // create a large enough buffer
    }

    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nUsers);  // count (broadcast) network uses of users
    double *Ecoop = calloc(sizeof(double), nUsers);  // record cooperation energy consumption of each user

    int decoded = 0;
    int count = 0;
    while (decoded != nUsers) {
        count++;
        pkt = snc_generate_packet(sc);

        // Phase I: broadcast to n users
        for (i=0; i<nUsers; i++) {
            if (snc_decoder_finished(decoders[i])) 
                continue;  // skip completed users

            nuse[i] += 1;
            if (rand() % 10000 >= delta[i] * 10000) {
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
            if (pD2D[k] != 0 && rand() % 10000 < pD2D[k]*10000) {
                // Cooperation among each pair
                for (int n=0; n<nD2Dmax[k]; n++) {
                    // U1 send to U2
                    if (n<nD2D[k*2] && !snc_decoder_finished(decoders[k*2+1])) {
                        pkt = snc_recode_packet(buffers[k*2], sched_t);
                        if (pkt == NULL)
                            continue;
                        Ecoop[k*2] += alpha;    // counting energy consumption for the cooperation packet
                        if (rand() % 10000 >= epsil[k]*10000) {
                            struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                            snc_process_packet(decoders[k*2+1], d2d_pkt);
                            snc_free_packet(d2d_pkt);
                        }
                        snc_free_packet(pkt);
                    }

                    // U2 send to U1
                    if (n<nD2D[k*2+1] && !snc_decoder_finished(decoders[k*2])) {
                        pkt = snc_recode_packet(buffers[k*2+1], sched_t);
                        if (pkt == NULL)
                            continue;
                        Ecoop[k*2+1] += alpha;    // counting energy consumption for the cooperation packet
                        if (rand() % 10000 >= epsil[k]*10000) {
                            struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                            snc_process_packet(decoders[k*2], d2d_pkt);
                            snc_free_packet(d2d_pkt);
                        }
                        snc_free_packet(pkt);
                    }
                }
            }
        }

        // check decoder status
        decoded = 0;
        for (i=0; i<nUsers; i++) {
            if (snc_decoder_finished(decoders[i]))
                decoded += 1;
        }
    }

    for (i=0; i<nUsers; i++) {
        struct snc_context *dsc = snc_get_enc_context(decoders[i]);
        unsigned char *rec_buf = snc_recover_data(dsc);
        if (memcmp(buf, rec_buf, sp.datasize) != 0)
            fprintf(stderr, "recovered is NOT identical to original.\n");
        printf("user: %d overhead: %.6f ops: %.6f network-uses: %d\n", i, snc_decode_overhead(decoders[i]), snc_decode_cost(decoders[i]), nuse[i]);
    }
    double totalEnergy = 0;
    for (i=0; i<npairs; i++) {
        printf("Pair-%d: ", i);
        for (j=0; j<2; j++)
            printf(" delta[%d]: %.3f ", j+1, delta[i*2+j]);
        printf(" pD2D: %.3f ", pD2D[i]);
        for (j=0; j<2; j++)
            printf(" nD2D[%d]: %d ", j+1, nD2D[i*2+j]);
        printf(" epsil[%d]: %.3f ", j+1, epsil[i]);
        for (j=0; j<2; j++)
            printf(" nuse[%d]: %d ", j+1, nuse[i*2+j]);
        printf(" CoopEnergy: %.6f ", Ecoop[i*2]+Ecoop[i*2+1]);
        totalEnergy += Ecoop[i*2]+Ecoop[i*2+1];
        printf("\n");
    }
    printf("snum: %d pktsize: %d GF_POWER: %d nBroadcast: %d TotalCoopEnergy: %.6f \n", sp.datasize/sp.size_p, sp.size_p, sp.gfpower, count, totalEnergy);
    snc_free_enc_context(sc);
    for (i=0; i<nUsers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}

