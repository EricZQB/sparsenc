#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "Simulate network coded multicasting with D2D\n\
                \n\
                \n\
                          ---> D1\n\
                S -----------> D2\n\
                          ---> ..\n\
                          ---> Dn\n\
                \n\
                usage: ./programName code_t decoder_t pktsize snum cnum size_b size_g delta epsil n_users D2Dinterval\n";

int main(int argc, char *argv[])
{
    if (argc != 12) {
        printf("%s\n", usage);
        exit(1);
    }
    int sched_t = MLPI_SCHED;

    struct snc_parameters sp;

    if (strcmp(argv[1], "RAND") == 0)
        sp.type = RAND_SNC;
    else if (strcmp(argv[1], "BAND") == 0)
        sp.type = BAND_SNC;
    else if (strcmp(argv[1], "WINDWRAP") == 0)
        sp.type = WINDWRAP_SNC;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    int decoder_t;
    if (strcmp(argv[2], "GG") == 0)
        decoder_t = GG_DECODER;
    else if (strcmp(argv[2], "OA") == 0)
        decoder_t = OA_DECODER;
    else if (strcmp(argv[2], "BD") == 0)
        decoder_t = BD_DECODER;
    else if (strcmp(argv[2], "CBD") == 0)
        decoder_t = CBD_DECODER;
    else if (strcmp(argv[2], "PP") == 0)
        decoder_t = PP_DECODER;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    sp.datasize = atoi(argv[3]) * atoi(argv[4]);
    sp.size_p   = atoi(argv[3]);
    sp.size_c   = atoi(argv[5]);
    sp.size_b   = atoi(argv[6]);
    sp.size_g   = atoi(argv[7]);
    sp.bpc      = 0;
    sp.bnc      = 0;
    sp.sys      = 0;
    sp.seed     = -1;

    double delta = atof(argv[8]);
    double epsil = atof(argv[9]);

    int nusers = atoi(argv[10]);
    int D2Dintval = atoi(argv[11]);

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
            if (rand() % 100 >= delta * 100) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                struct snc_packet *buf_pkt = snc_duplicate_packet(pkt, &sp);
                snc_buffer_packet(buffers[i], buf_pkt);              // save a copy in the buffer, for the purpose of D2D comm.
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
        
        // Phase II: D2D communication
        // Each user send a packet to other users
        // Do D2D communication every 10 time slots
        if (D2Dintval != 0 && count % D2Dintval == 0) {
            for (i=0; i<nusers; i++) {
                pkt = snc_recode_packet(buffers[i], sched_t);
                if (pkt == NULL)
                    continue;
                
                for (j=0; j<nusers; j++) {
                    if (j==i || snc_decoder_finished(decoders[j]))
                        continue;  // skip himself and decoded
                    // D2D from i to j
                    if (rand() % 100 >= epsil * 100) {
                        struct snc_packet *d2d_pkt = snc_duplicate_packet(pkt, &sp);
                        snc_process_packet(decoders[j], d2d_pkt);
                        snc_free_packet(d2d_pkt);
                    }
                }
                snc_free_packet(pkt);
            }
        }

        // check decoder status
        decoded = 0;
        for (i=0; i<nusers; i++) {
            if (snc_decoder_finished(decoders[i]))
                decoded += 1;
        }
    }

    printf("nusers: %d snum: %d pktsize: %d delta: %.3f epsil: %.3f\n", nusers, sp.datasize/sp.size_p, sp.size_p, delta, epsil);
    for (i=0; i<nusers; i++) {
        struct snc_context *dsc = snc_get_enc_context(decoders[i]);
        unsigned char *rec_buf = snc_recover_data(dsc);
        if (memcmp(buf, rec_buf, sp.datasize) != 0)
            fprintf(stderr, "recovered is NOT identical to original.\n");
        printf("user: %d overhead: %.6f ops: %.6f network-uses: %d\n", i, snc_decode_overhead(decoders[i]), snc_decode_cost(decoders[i]), nuse[i]);
    }
    snc_free_enc_context(sc);
    for (i=0; i<nusers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}
