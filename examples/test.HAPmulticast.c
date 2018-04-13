#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "Simulate systematic network coded multicasting over an HAP relay\n\
                \n\
                \n\
                          ---> D1\n\
                S --- HAP ---> D2\n\
                          ---> ..\n\
                          ---> Dn\n\
                \n\
                usage: ./programName snum size_p bufsize delta epsil n_users\n";

int main(int argc, char *argv[])
{
    if (argc != 7) {
        printf("%s\n", usage);
        exit(1);
    }
    int decoder_t = CBD_DECODER;
    int sched_t = RAND_SCHED_SYS;

    struct snc_parameters sp;
    sp.type = BAND_SNC;
    sp.datasize = atoi(argv[1]) * atoi(argv[2]);
    sp.size_p   = atoi(argv[2]);
    sp.size_c   = 0;
    sp.size_b   = atoi(argv[1]);
    sp.size_g   = atoi(argv[1]);
    sp.bpc      = 0;
    sp.bnc      = 0;
    sp.sys      = 1;
    sp.seed     = -1;

    int bufsize = atoi(argv[3]);
    double delta = atof(argv[4]);
    double epsil = atof(argv[5]);

    int nusers = atoi(argv[6]);

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

    /* Create NC buffer at HAP*/
    struct snc_buffer *buffer ;
    if ((buffer = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer.\n");
        exit(1);
    }

    /* Create decoders */
    sp.seed = (snc_get_parameters(sc))->seed;
    struct snc_decoder **decoders = malloc(sizeof(struct snc_decoder *) * nusers);
    for (i=0; i<nusers; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
    }
    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nusers);  // count network uses of each user

    int decoded = 0;
    while (decoded != nusers) {
        pkt = snc_generate_packet(sc);
        
        if (rand() % 100 >= delta * 100) {
            snc_buffer_packet(buffer, pkt);
        } else {
            snc_free_packet(pkt);
        }

        // generate a recoded packet and broadcast
        pkt = snc_alloc_empty_packet(&sp);
        snc_recode_packet_im(buffer, pkt, sched_t);
        for (i=0; i<nusers; i++) {
            if (!snc_decoder_finished(decoders[i])) {
                struct snc_packet *dup_pkt = snc_duplicate_packet(pkt, &sp);
                nuse[i] += 1;
                if (rand() % 100 >= epsil * 100) {
                    snc_process_packet(decoders[i], dup_pkt);
                }
                snc_free_packet(dup_pkt);

                if (snc_decoder_finished(decoders[i]))
                    decoded += 1;
            }
        }
        snc_free_packet(pkt);
    }

    printf("bufsize: %d nusers: %d snum: %d pktsize: %d delta: %.3f epsil: %.3f\n", bufsize, nusers, sp.datasize/sp.size_p, sp.size_p, delta, epsil);
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
