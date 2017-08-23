#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "Simulate lossy butterfly networks:\n\
                    A-----------------(R1) \n\
                   / \\               /      \n\
                  /   \\             /      \n\
                 S     C------------D      \n\
                  \\   /             \\     \n\
                   \\ /               \\    \n\
                    B-----------------(R2) \n\
                We only need to simulate for R1. Transmissions on B->R2 and D->R2 are omitted.\n\
                \n\
                usage: ./programName code_t dec_t sched_t datasize size_p\n\
                                     size_c size_b size_g bpc bnc sys bufsize\n\
                                     pe1 pe2 ...\n\
                code_t  - code type: RAND, BAND, WINDWRAP\n\
                dec_t   - decoder type: GG, OA, BD, CBD, PP\n\
                sched_t - scheduling type: TRIV, RAND, MLPI, NURAND\n\
                datasize - Number of bytes\n\
                size_p   - Packet size in bytes\n\
                size_c   - Nnumebr of check packets\n\
                size_b   - Subgeneration distance\n\
                size_g   - Subgeneration size\n\
                bpc      - Use binary precode (0 or 1)\n\
                bnc      - Use binary network code (0 or 1)\n\
                sys      - Systematic code (0 or 1)\n\
                bufsize  - buffer size of each intermediate node\n\
                pe_i     - erasure probabilities of SA, SB, AC, BC, AR, CD, DR, respectively.\n\
                           If only one pe_i is given, all the links have the same erasure rate.\n";
int main(int argc, char *argv[])
{
    if (argc != 14 && argc != 20) {
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
    else if (strcmp(argv[3], "MLPI") == 0)
        sched_t = MLPI_SCHED;
    else if (strcmp(argv[3], "NURAND") == 0)
        sched_t = NURAND_SCHED;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    sp.datasize = atoi(argv[4]);
    sp.size_p   = atoi(argv[5]);
    sp.size_c   = atoi(argv[6]);
    sp.size_b   = atoi(argv[7]);
    sp.size_g   = atoi(argv[8]);
    sp.bpc      = atoi(argv[9]);
    sp.bnc      = atoi(argv[10]);
    sp.sys      = atoi(argv[11]);
    sp.seed     = -1;
    int bufsize = atoi(argv[12]);
    double Esa, Esb, Eac, Ebc, Ear, Ecd, Edr;
    if (argc == 14)
        Esa = Esb = Eac = Ebc = Ear = Ecd = Edr = atof(argv[13]);
    else {
        Esa = atof(argv[13]);
        Esb = atof(argv[14]);
        Eac = atof(argv[15]);
        Ebc = atof(argv[16]);
        Ear = atof(argv[17]);
        Ecd = atof(argv[18]);
        Edr = atof(argv[19]);
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
    struct snc_buffer *buf_A;
    if ((buf_A = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer for node A.\n");
        exit(1);
    }
    struct snc_buffer *buf_B;
    if ((buf_B = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer for node B.\n");
        exit(1);
    }
    struct snc_buffer *buf_C;
    if ((buf_C = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer for node C.\n");
        exit(1);
    }
    struct snc_buffer *buf_D;
    if ((buf_D = snc_create_buffer(snc_get_parameters(sc), bufsize)) == NULL) {
        fprintf(stderr, "Cannot create snc buffer for node D.\n");
        exit(1);
    }

    /* Create decoder */
    sp.seed = (snc_get_parameters(sc))->seed;
    struct snc_decoder *decoder = snc_create_decoder(&sp, decoder_t);
    clock_t start, stop, dtime = 0;
    clock_t recv_start;
    int nuse = 0;
    int first_recv = 1;
    while (!snc_decoder_finished(decoder)) {
        nuse++;
        // S-->A
        struct snc_packet *pktSA = snc_generate_packet(sc);
        if (rand() % 100 >= Esa * 100)
            snc_buffer_packet(buf_A, pktSA);
        else
            snc_free_packet(pktSA);
        // S-->B
        struct snc_packet *pktSB = snc_generate_packet(sc);
        if (rand() % 100 >= Esb * 100)
            snc_buffer_packet(buf_B, pktSB);
        else
            snc_free_packet(pktSB);
        // A-->C
        struct snc_packet *pktAC = snc_alloc_empty_packet(&sp);
        if (snc_recode_packet_im(buf_A, pktAC, sched_t) != -1 && rand()%100 >= Eac*100) {
            snc_buffer_packet(buf_C, pktAC);
        } else {
            snc_free_packet(pktAC);
        }
        // B-->C
        struct snc_packet *pktBC = snc_alloc_empty_packet(&sp);
        if (snc_recode_packet_im(buf_B, pktBC, sched_t) != -1 && rand()%100 >= Ebc*100) {
            snc_buffer_packet(buf_C, pktBC);
        } else {
            snc_free_packet(pktBC);
        }
        // C-->D
        struct snc_packet *pktCD = snc_alloc_empty_packet(&sp);
        if (snc_recode_packet_im(buf_C, pktCD, sched_t) != -1 && rand()%100 >= Ecd*100) {
            snc_buffer_packet(buf_D, pktCD);
        } else {
            snc_free_packet(pktCD);
        }

        // A-->R, R performs on-the-fly decoding
        struct snc_packet *pktAR = snc_alloc_empty_packet(&sp);
        if (snc_recode_packet_im(buf_A, pktAR, sched_t) != -1 && rand()%100 >= Ear*100) {
            /* Measure decoding time */
            start = clock();
            if (first_recv) {
                recv_start = clock();
                first_recv = 0;
            }
            snc_process_packet(decoder, pktAR);
            stop = clock();
            dtime += stop - start;
        }
        snc_free_packet(pktAR);

        // D-->R, R performs on-the-fly decoding
        struct snc_packet *pktDR = snc_alloc_empty_packet(&sp);
        if (snc_recode_packet_im(buf_D, pktDR, sched_t) != -1 && rand()%100 >= Edr*100) {
            /* Measure decoding time */
            start = clock();
            snc_process_packet(decoder, pktDR);
            stop = clock();
            dtime += stop - start;
        }
        snc_free_packet(pktDR);
    }
    clock_t decode_delay = clock() - recv_start;
    printf("dec-time: %.6f dec-delay: %.6f bufsize: %d Esa: %.3f Esb: %.3f Eac: %.3f Ebc: %.3f Ear: %.3f Ecd: %.3f Edr: %.3f nuses: %d\n", ((double) dtime)/CLOCKS_PER_SEC, ((double) decode_delay)/CLOCKS_PER_SEC, bufsize, Esa, Esb, Eac, Ebc, Ear, Ecd, Edr, nuse);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");
    print_code_summary(dsc, snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_buffer(buf_A);
    snc_free_buffer(buf_B);
    snc_free_buffer(buf_C);
    snc_free_buffer(buf_D);
    snc_free_decoder(decoder);
    return 0;
}
