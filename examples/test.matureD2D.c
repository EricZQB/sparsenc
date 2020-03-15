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
                usage: ./programName nUsers headUserIdx B2U_loss_probabilities H2U_loss_probabilities\n";

int main(int argc, char *argv[])
{
    if (argc < 4) {
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
    int headIdx = atoi(argv[2]);                    // designated head user (count from 1 as it is generated from Matlab)
    double *delta = calloc(nUsers, sizeof(double)); // broadcast link loss probability
    double *epsil = calloc(nUsers, sizeof(double)); // loss probabilities from the head user to each of the others
    int i, j;
    for (i=0; i<nUsers; i++) {
        delta[i] = atof(argv[3+i]);
        epsil[i] = atof(argv[3+nUsers+i]);
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

    // Create decoders and buffer for each user of each pair
    struct snc_decoder **decoders =  malloc(sizeof(struct snc_decoder*) * nUsers);
    for (i=0; i<nUsers; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
    }

    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nUsers);  // count network uses of users
    double *Ecoop = calloc(sizeof(double), nUsers);  // record cooperation energy consumption of each user

    int count1 = 0;

    // Phase I, broadcast to all users
    while ( !snc_decoder_finished(decoders[headIdx-1]) ) {
        count1++;
        pkt = snc_generate_packet(sc);

        // Phase I: broadcast to n users
        for (i=0; i<nUsers; i++) {
            if (snc_decoder_finished(decoders[i])) 
                continue;  // skip completed users

            nuse[i] += 1;
            if (rand() % 10000 >= delta[i] * 10000) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
    }
    printf("Head user is complete after %d broadcasts \n", count1);
    // check decoder status
    int decoded = 0;
    for (i=0; i<nUsers; i++) {
        if (snc_decoder_finished(decoders[i]))
            decoded += 1;
    }
    int count2 = 0;
    // Phase II, broadcast from the completed user (head) to others
    // Since it's the same file, we can just re-use sc to simulate generating coded packets from the head user
    while (decoded != nUsers) {
        count2++;
        pkt = snc_generate_packet(sc);

        // Phase I: broadcast to n users
        for (i=0; i<nUsers; i++) {
            if (i==headIdx-1 || snc_decoder_finished(decoders[i])) 
                continue;  // skip completed users

            nuse[i] += 1;
            if (rand() % 10000 >= epsil[i] * 10000) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);

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
    printf("snum: %d pktsize: %d GF_POWER: %d nBsBroadcast: %d nHeadBroadcast: %d \n", sp.datasize/sp.size_p, sp.size_p, sp.gfpower, count1, count2);
    snc_free_enc_context(sc);
    for (i=0; i<nUsers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}

