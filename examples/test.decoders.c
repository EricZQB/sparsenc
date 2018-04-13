#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

char usage[] = "usage: ./sncDecoder code_t dec_t datasize size_p size_c size_b size_g bpc bnc sys\n\
                       code_t   - RAND, BAND, WINDWRAP\n\
                       dec_t    - GG, OA, BD, CBD, PP\n\
                       datasize - Number of bytes\n\
                       size_p   - Packet size in bytes\n\
                       size_c   - Nnumebr of check packets\n\
                       size_b   - Subgeneration distance\n\
                       size_g   - Subgeneration size\n\
                       bpc      - Use binary precode (0 or 1)\n\
                       bnc      - Use binary network code (0 or 1)\n\
                       sys      - Systematic code (0 or 1)\n";
int main(int argc, char *argv[])
{
    if (argc != 11) {
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
    else if (strcmp(argv[1], "BATS") == 0)
        sp.type = BATS_SNC;
    else {
        printf("%s\n", usage);
        exit(1);
    }

    int decoder_type;
    if (strcmp(argv[2], "GG") == 0)
        decoder_type = GG_DECODER;
    else if (strcmp(argv[2], "OA") == 0)
        decoder_type = OA_DECODER;
    else if (strcmp(argv[2], "BD") == 0)
        decoder_type = BD_DECODER;
    else if (strcmp(argv[2], "CBD") == 0)
        decoder_type = CBD_DECODER;
    else if (strcmp(argv[2], "PP") == 0)
        decoder_type = PP_DECODER;
    else {
        printf("%s\n", usage);
        exit(1);
    }
    sp.datasize = atoi(argv[3]);
    sp.size_p   = atoi(argv[4]);
    sp.size_c   = atoi(argv[5]);
    sp.size_b   = atoi(argv[6]);
    sp.size_g   = atoi(argv[7]);
    sp.bpc      = atoi(argv[8]);
    sp.bnc      = atoi(argv[9]);
    sp.sys      = atoi(argv[10]);
    sp.seed     = -1;  // Initialize seed as -1

    char *ur = getenv("SNC_NONUNIFORM_RAND");
    if ( ur != NULL && atoi(ur) == 1) {
        if (sp.type != BAND_SNC || sp.size_b != 1) {
            printf("Non-Uniform Random Scheduling can only be used for BAND code with size_b=1.\n");
            exit(1);
        }
    }

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
    // Make decoder stop in the middle of decoding.
    // Test saving/restoring decoder context to/from file.
    /*
    int count = 0;
    while (snc_decoder_finished(decoder) != 1) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        //Measure decoding time
        start = clock();
        snc_process_packet(decoder, pkt);
        stop = clock();
        dtime += stop - start;
        count++;
        if (count > 2200) {
            printf("Save decoder context into file\n");
            snc_save_decoder_context(decoder, "CBDdecoder.part");
            break;
        }
    }
    snc_free_decoder(decoder);

    decoder = snc_restore_decoder("CBDdecoder.part");
    */
    clock_t send_start = clock();
    clock_t decode_delay = 0;
    while (snc_decoder_finished(decoder) != 1) {
        struct snc_packet *pkt = snc_generate_packet(sc);
        /* Measure decoding time */
        start = clock();
        snc_process_packet(decoder, pkt);
        snc_free_packet(pkt);
        stop = clock();
        dtime += (stop - start);
    }
    decode_delay = clock() - send_start;
    //printf("clocks: %d CLOCKS_PER_SEC: %d \n", dtime, CLOCKS_PER_SEC);

    printf("dec-time: %.6f decode-delay: %.6f ", (double) dtime/CLOCKS_PER_SEC, (double) decode_delay/CLOCKS_PER_SEC);

    struct snc_context *dsc = snc_get_enc_context(decoder);
    unsigned char *rec_buf = snc_recover_data(dsc);
    if (memcmp(buf, rec_buf, sp.datasize) != 0)
        fprintf(stderr, "recovered is NOT identical to original.\n");

    print_code_summary(dsc, snc_decode_overhead(decoder), snc_decode_cost(decoder));

    snc_free_enc_context(sc);
    snc_free_decoder(decoder);
    return 0;
}
