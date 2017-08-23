#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "usage: ./sncRLNC pktsize pktnum bnc\n\
                       pktsize  - Packet size in bytes\n\
                       pktnum   - Numebr of packets\n\
                       bnc      - Use binary network code (0 or 1), default 0\n";

int main(int argc, char *argv[])
{
    if (argc != 3 && argc != 4) {
        printf("%s\n", usage);
        exit(1);
    }
    struct snc_parameters sp;
    sp.type = BAND_SNC;

    int decoder_type = CBD_DECODER;
    sp.size_p   = atoi(argv[1]);
    sp.datasize = atoi(argv[2]) * atoi(argv[1]);
    sp.size_c   = 0;
    sp.size_b   = atoi(argv[2]);
    sp.size_g   = atoi(argv[2]);
    sp.bpc      = 0;
    sp.sys      = 0;
    sp.seed     = -1;  // Initialize seed as -1
    if (argc == 4 && atoi(argv[3]) == 1)
        sp.bnc = 1;
    else
        sp.bnc = 0;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    unsigned char *buf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    if ((sc = snc_create_enc_context(buf, &sp)) == NULL) {
        fprintf(stderr, "Cannot create File Context.\n");
        return 1;
    }

    sp.seed = (snc_get_parameters(sc))->seed;
    struct snc_decoder *decoder = snc_create_decoder(&sp, decoder_type);
    if (decoder == NULL)
        exit(1);
    clock_t start, stop, dtime = 0;

    while (snc_decoder_finished(decoder) != 1) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        /* Measure decoding time */
        start = clock();
        snc_process_packet(decoder, pkt);
        snc_free_packet(pkt);
        stop = clock();
        dtime += (stop - start);
    }
    //printf("clocks: %d CLOCKS_PER_SEC: %d \n", dtime, CLOCKS_PER_SEC);

    printf("dec-time: %.6f ", (double) dtime/CLOCKS_PER_SEC);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");

    printf("pktnum: %d pktsize: %d overhead: %.3f complexity: %.3f (operation per information symbol)\n", atoi(argv[2]), atoi(argv[1]), snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_decoder(decoder);
    return 0;
}
