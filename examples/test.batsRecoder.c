#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#include "sparsenc.h"

char usage[] = "Simulate fiex-degree BATS codes over n-hop lossy line networks\n\
                \n\
                S ---(1-pe_1)---> V1 ---(1-pe_2)---> ... ---(1-pe_n)---> D\n\
                \n\
                usage: ./programName density pktsize snum cnum bts bufsize nhop pe1 pe2 R1 R2...\n\
                density  - maximum allowed degree of each coded packet\n\
                pktsize  - packet size (in bytes)\n\
                snum     - number of source packets\n\
                cnum     - number of parity-check packets of precode\n\
                bts      - batch transmission size: # of packets sent for each batch\n\
                bufsize  - buffer size at intermediate nodes\n\
                nhop     - number of hops (integer)\n\
                pe_i     - erasure probabilities of each hop. The number of pe_i's is either equal to nhop,\n\
                           or only one which corresponds to the homogeneous case (all hops have the same erausre rate)\n\
                R_i      - Number of packets sent on each hop in each time slot, usage is similar as pe_i.\n";
int main(int argc, char *argv[])
{
    if (argc < 10 || (argc != 10 && argc != 8 + 2 * atoi(argv[7]))) {
        printf("%s\n", usage);
        exit(1);
    }

    int i, j;
    int gsize = atoi(argv[1]);
    int pktsize = atoi(argv[2]);
    int snum = atoi(argv[3]);
    int cnum = atoi(argv[4]);
    int bsize = atoi(argv[5]);
    int bufsize = atoi(argv[6]);
    int numhop = atoi(argv[7]);
    int datasize = snum * pktsize;
    
    double *pe = calloc(numhop, sizeof(double));
    int *rate = calloc(numhop, sizeof(int));
    for (i=0; i<numhop; i++) {
        if (argc == 10) {
            pe[i] = atof(argv[8]);
            rate[i] = atoi(argv[9]);
        } else {
            pe[i] = atof(argv[8+i]);
            rate[i] = atoi(argv[8+numhop+i]);
        }
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    // generate random source data
    unsigned char *databuf = malloc(datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, databuf, datasize);
    close(rnd);

    struct snc_parameters sp;

    sp.datasize = datasize;
    sp.size_p   = pktsize;
    sp.size_c   = cnum;
    sp.size_b   = bsize;
    sp.size_g   = gsize;
    sp.bpc      = 1;
    sp.bnc      = 0;
    sp.sys      = 0;
    sp.seed     = -1;
    sp.type     = BATS_SNC;


    struct snc_context *sc;
    /* Create GNC encoding context */
    if ((sc = snc_create_enc_context(databuf, &sp)) == NULL) {
        fprintf(stderr, "Cannot create snc_context.\n");
        exit(1);
    }

    /* Create recoder buffers */
    // n-hop network has (n-1) intermediate nodes, and therefore has (n-1) recoders
    struct snc_buffer_bats **buffer = malloc(sizeof(struct snc_buffer_bats*) * (numhop-1));
    for (i=0; i<numhop-1; i++) {
        if ((buffer[i] = snc_create_buffer_bats(snc_get_parameters(sc), bufsize)) == NULL) {
            fprintf(stderr, "Cannot create snc buffer.\n");
            exit(1);
        }
    }

    /* Create decoder */
    int decoder_t = OA_DECODER;
    sp.seed = (snc_get_parameters(sc))->seed;
    struct snc_decoder *decoder = snc_create_decoder(&sp, decoder_t);
    clock_t start, stop, dtime = 0;
    clock_t recv_start;
    struct snc_packet *pkt;    // pointer of coded packet
    int nuse = 0;  // count network uses
    int first_recv = 1;
    while (!snc_decoder_finished(decoder)) {
        nuse++;
        for (i=0; i<numhop; i++) {
            for (j=0; j<rate[i]; j++) {
                if (i == 0) {
                    pkt = snc_generate_packet(sc);  // coded packet generated at the source node
                } else {
                    pkt = snc_alloc_empty_packet(&sp);
                    if (snc_recode_packet_bats_im(buffer[i-1], pkt) == -1)
                        continue;
                }
                if (rand() % 100 >= pe[i] * 100) {
                    if (i < numhop-1) {
                        snc_buffer_packet_bats(buffer[i], pkt);   // intermediate ndoes buffer packets
                    } else {
                        /* Measure decoding time */
                        start = clock();
                        if (first_recv) {
                            recv_start = clock();
                            first_recv = 0;
                        }
                        snc_process_packet(decoder, pkt);
                        stop = clock();
                        dtime += stop - start;
                        snc_free_packet(pkt);
                    }
                } else {
                    snc_free_packet(pkt);
                }
            }
        }
    }
    clock_t decode_delay = clock() - recv_start;
    printf("dec-time: %.6f dec-delay: %.6f bufsize: %d numhop: %d ", ((double) dtime)/CLOCKS_PER_SEC, (double) decode_delay/CLOCKS_PER_SEC, bufsize, numhop);
    for (i=0; i<numhop; i++) {
        printf("rate[%i]: %d ", i, rate[i]);
        printf("pe[%i]: %.3f ", i, pe[i]);
    }
    printf(" nuses: %d\n", nuse);


    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(databuf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");
    print_code_summary(dsc, snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    for (i=0; i<numhop-1; i++) 
        snc_free_buffer_bats(buffer[i]);
    snc_free_decoder(decoder);
    return 0;
}
