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
                usage: GF_SIZE=q ./programName pktsize snum delta1 delta2 epsil1 epsil2 pD2D nD2D_1 nD2D_2\n";

int main(int argc, char *argv[])
{
    if (argc != 10) {
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
    sp.sys      = 0;
    sp.seed     = -1;

	double delta[2] = {atof(argv[3]), atof(argv[4])};
    double epsil[2] = {atof(argv[5]), atof(argv[6])};

    int nusers = 2;
	
    double pD2D = atof(argv[7]);
    int nD2D_1 = atoi(argv[8]);
    int nD2D_2 = atoi(argv[9]);
    int nD2D = nD2D_1 > nD2D_2 ? nD2D_1 : nD2D_2;
	
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

    // Create decoders and buffer for each user
    struct snc_decoder **decoders =  malloc(sizeof(struct snc_decoder*) * nusers);
    struct snc_buffer **buffers = malloc(sizeof(struct snc_buffer*) * nusers);
    for (i=0; i<nusers; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
        buffers[i] = snc_create_buffer(&sp, sp.size_g);  // create a large enough buffer
    }

    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nusers);  // count network uses of each user

    int decoded = 0;
    int count = 0;
    while (decoded != nusers) {
        count++;
        pkt = snc_generate_packet(sc);

        // Phase I: broadcast to n users
        for (i=0; i<nusers; i++) {
            if (snc_decoder_finished(decoders[i])) 
                continue;

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
        if (pD2D != 0 && rand() % 100 < pD2D*100) {
            for (int n=0; n<nD2D; n++) {
                // U1 send to U2
                if (n<nD2D_1 && !snc_decoder_finished(decoders[1])) {
                    pkt = snc_recode_packet(buffers[0], sched_t);
                    if (pkt == NULL)
                        continue;
                    if (rand() % 100 >= epsil[0]) {
                        struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                        snc_process_packet(decoders[1], d2d_pkt);
                        snc_free_packet(d2d_pkt);
                    }
                    snc_free_packet(pkt);
                }

                // U2 send to U1
                if (n<nD2D_2 && !snc_decoder_finished(decoders[0])) {
                    pkt = snc_recode_packet(buffers[1], sched_t);
                    if (pkt == NULL)
                        continue;
                    if (rand() % 100 >= epsil[1]) {
                        struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                        snc_process_packet(decoders[0], d2d_pkt);
                        snc_free_packet(d2d_pkt);
                    }
                    snc_free_packet(pkt);
                }
            }
	    }

        // check decoder status
        decoded = 0;
        for (i=0; i<nusers; i++) {
            if (snc_decoder_finished(decoders[i]))
                decoded += 1;
        }
    }

    for (i=0; i<nusers; i++) {
        struct snc_context *dsc = snc_get_enc_context(decoders[i]);
        unsigned char *rec_buf = snc_recover_data(dsc);
        if (memcmp(buf, rec_buf, sp.datasize) != 0)
            fprintf(stderr, "recovered is NOT identical to original.\n");
        printf("user: %d overhead: %.6f ops: %.6f network-uses: %d\n", i, snc_decode_overhead(decoders[i]), snc_decode_cost(decoders[i]), nuse[i]);
    }
    int GF_size = snc_get_GF_size(&sp);
    printf("snum: %d pktsize: %d delta1: %.3f delta2: %.3f epsil1: %.3f epsil2: %.3f nuse1: %d nuse2: %d GF_SIZE: %d pD2D: %.3f nD2D_1: %d nD2D_2: %d \n", sp.datasize/sp.size_p, sp.size_p, delta[0], delta[1], epsil[0], epsil[1], nuse[0], nuse[1], GF_size, pD2D, nD2D_1, nD2D_2);
    snc_free_enc_context(sc);
    for (i=0; i<nusers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}

