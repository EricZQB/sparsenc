#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "Simulate n-hop lossy line networks\n\
                \n\
                S ---R1(1-pe_1)---> V1 ---R2(1-pe_2)---> ... ---Rn(1-pe_n)---> D\n\
                \n\
                usage: ./programName code_t dec_t sched_t datasize size_p size_c\n\
                                     size_b size_g bpc gfpower sys bufsize\n\
                                     nhop R1 R2 ... pe1 pe2 ...\n\
                code_t  - code type: RAND, BAND, WINDWRAP\n\
                dec_t   - decoder type: GG, OA, BD, CBD\n\
                sched_t - scheduling type: TRIV, RAND, RANDSYS, MLPI, MLPISYS, NURAND\n\
                datasize - bytes of data to send\n\
                size_p   - packet size (in bytes)\n\
                size_c   - number of parity-check packets of precode\n\
                size_b   - base subgeneration size\n\
                size_g   - subgeneration size (after adding overlap)\n\
                bpc      - Use binary precode (0 or 1)\n\
                gfpower  - power of Galois field (1,2,...,8)\n\
                sys      - Systematic code (0 or 1)\n\
                bufsize  - buffer size of each intermediate node\n\
                nhop     - number of hops (integer)\n\
                R_i      - packet sending rate (integer) of each hop. The number of R_i's is either equal to nhop,\n\
                           or 1 which corresponds to that all the hops have equal packet sending rate.\n\
                pe_i     - erasure probabilities of each hop. The number of pe_i's is either equal to nhop,\n\
                           or 1 which corresponds to the homogeneous case (all hops have the same erausre rate)\n";
int main(int argc, char *argv[])
{
    if (argc < 15 || (argc != 16 && argc != 15 + atoi(argv[13]) && argc != 14+atoi(argv[13])*2)) {
        printf("%s\n", usage);
        exit(1);
    }
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

    int sched_t;
    if (strcmp(argv[3], "TRIV") == 0)
        sched_t = TRIV_SCHED;
    else if (strcmp(argv[3], "RAND") == 0)
        sched_t = RAND_SCHED;
    else if (strcmp(argv[3], "RANDSYS") == 0)
        sched_t = RAND_SCHED_SYS;
    else if (strcmp(argv[3], "MLPI") == 0)
        sched_t = MLPI_SCHED;
    else if (strcmp(argv[3], "MLPISYS") == 0)
        sched_t = MLPI_SCHED_SYS;
    else if (strcmp(argv[3], "NURAND") == 0)
        sched_t = NURAND_SCHED;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    sp.datasize = atoi(argv[4]);
    sp.size_p   = atof(argv[5]);
    sp.size_c   = atoi(argv[6]);
    sp.size_b   = atoi(argv[7]);
    sp.size_g   = atoi(argv[8]);
    sp.bpc      = atoi(argv[9]);
    sp.gfpower  = atoi(argv[10]);
    sp.sys      = atoi(argv[11]);
    sp.seed     = -1;
    int bufsize = atoi(argv[12]);
    int numhop  = atoi(argv[13]);    // Number of hops of the line network
    int *rate   = malloc(sizeof(int) * numhop);
    double *pe  = malloc(sizeof(double) * numhop);
    int i, j;
    for (i=0; i<numhop; i++) {
        if (argc == 16) {
            rate[i] = atoi(argv[14]);
            pe[i] = atof(argv[15]);
        } else if (argc == 14 + atoi(argv[13]) + 1 && atof(argv[15]) >= 1) {
            rate[i] = atoi(argv[14+i]);
            pe[i] = atof(argv[14+numhop]);
        } else if (argc == 14 + atoi(argv[13]) + 1 && atof(argv[15]) < 1) {
            rate[i] = atoi(argv[14]);
            pe[i] = atof(argv[15+i]);        
        } else {
            rate[i] = atoi(argv[14+i]);
            pe[i]   = atof(argv[14+numhop+i]);
        }
    }

    char *ur = getenv("SNC_NONUNIFORM_RAND");
    if ( ur != NULL && atoi(ur) == 1) {
        if (sp.type != BAND_SNC || sp.size_b != 1) {
            printf("Non-Uniform Random Scheduling can only be used for BAND code with size_b=1.\n");
            exit(1);
        }
    }

    if (sched_t == NURAND_SCHED) {
        if (sp.type != BAND_SNC || sp.size_b != 1) {
            printf("Non-Uniform Random Scheduling can only be used for BAND code with size_b=1.\n");
            exit(1);
        }
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

    /* Create recoder buffers */
    // n-hop network has (n-1) intermediate nodes, and therefore has (n-1) recoders
    struct snc_buffer **buffer = malloc(sizeof(struct snc_buffer*) * (numhop-1));
    for (i=0; i<numhop-1; i++) {
        if ((buffer[i] = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
            fprintf(stderr, "Cannot create snc buffer.\n");
            exit(1);
        }
    }

    /* Create decoder */
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
                    if (snc_recode_packet_im(buffer[i-1], pkt, sched_t) == -1)
                        continue;
                }
                if (rand() % 10000000 >= pe[i] * 10000000) {
                    if (i < numhop-1) {
                        snc_buffer_packet(buffer[i], pkt);   // intermediate ndoes buffer packets
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
    int GF_size = snc_get_GF_power(&sp);
    printf("GF_SIZE: %d \n", (1<<GF_size));
    printf("dec-time: %.6f dec-delay: %.6f bufsize: %d numhop: %d ", ((double) dtime)/CLOCKS_PER_SEC, (double) decode_delay/CLOCKS_PER_SEC, bufsize, numhop);
    for (i=0; i<numhop; i++) {
        printf("rate[%i]: %d ", i, rate[i]);
        printf("pe[%i]: %.7f ", i, pe[i]);
    }
    printf(" nuses: %d\n", nuse);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");
    print_code_summary(dsc, snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    for (i=0; i<numhop-1; i++) 
        snc_free_buffer(buffer[i]);
    snc_free_decoder(decoder);
    return 0;
}
