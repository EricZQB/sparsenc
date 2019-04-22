#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

// Gilbert-Elliott model
#define CH_GOOD 0
#define CH_BAD  1
int channel_state[256] = {CH_GOOD};  // save channel states for up to 256 links

char usage[] = "Simulate n-hop lossy line networks\n\
                \n\
                S ---R1(1-pe_1)---> V1 ---R2(1-pe_2)---> ... ---Rn(1-pe_n)---> D\n\
                \n\
                usage: ./programName code_t dec_t sched_t datasize size_p size_c\n\
                                     size_b size_g bpc gfpower sys bufsize\n\
                                     nhop pg pb ag ab\n\
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
                pg       - probability of transition from Good to Bad state\n\
                pb       - probability of transition from Bad to Good state\n\
                ag       - Proability of successful transmission in the Good state\n\
                ab       - Proability of successful transmission in the Bad state\n";
int main(int argc, char *argv[])
{
    if (argc < 18 || (argc != 18 && argc != 14 + 4*atoi(argv[13]))) {
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
    int *rate = calloc(numhop, sizeof(int));

    int i, j;
    // (pg, pb, ag, ab) probabilities
    double *prob = calloc(numhop*4, sizeof(double));
    for (i=0; i<numhop; i++) {
        rate[i] = 1;
        if (argc == 18) {
            // Each hop has identical parameters (pg, pb, h)
            prob[i*4] = atof(argv[14]);
            prob[i*4+1] = atof(argv[15]);
            prob[i*4+2] = atof(argv[16]);
            prob[i*4+3] = atof(argv[17]);
        } else {
            prob[i*4] = atof(argv[14+i*4]);
            prob[i*4+1] = atof(argv[15+i*4]);
            prob[i*4+2] = atof(argv[16+i*4]);
            prob[i*4+3] = atof(argv[17+i*4]);
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
                // Gilbert packet loss model
                if (!Gilbert_Elliott_erasure(prob[4*i], prob[i*4+1], prob[i*4+2], prob[i*4+3], i)) {
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
        printf("pg[%i]: %.3f ", i, prob[i*4]);
        printf("pb[%i]: %.3f ", i, prob[i*4+1]);
        printf("ag[%i]: %.3f ", i, prob[i*4+2]);
        printf("ab[%i]: %.3f ", i, prob[i*4+3]);
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


// Determine whether packet loss occurs for link i
int Gilbert_Elliott_erasure(double pg, double pb, double ag, double ab, int i)
{
    // Determine channel state
    if (channel_state[i] == CH_GOOD) {
        if (rand() % 100 < pg * 100) {
            // transition from Good to Bad
            channel_state[i] = CH_BAD;
        }
    } else {
        if (rand() % 100 < pb * 100) {
            // transition from Bad to Good
            channel_state[i] = CH_GOOD;
        }
    }

    // Determine erasure depending on the channel state
    // Erasure occurs with probability 1-h in the Bad state, and with probability 0 in the Good state
    if ((channel_state[i] == CH_GOOD && rand() % 100 < (1-ag) * 100) || (channel_state[i] == CH_BAD && rand() % 100 < (1-ab) * 100))
        return 1;
    else
        return 0;
}
