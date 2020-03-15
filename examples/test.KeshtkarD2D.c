#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sparsenc.h"

// Simulate Keshtkarjahromi-TNET-2018 NCMI-batch scheme (algorithm 2)
char usage[] = "Simulate RLNC network coded D2D cooperative transmissions\n\
                \n\
                \n\
                    --------U1\n\
                   /        ||\n\
                S -         ||\n\
                   \\       ||\n\
                    --------U2\n\
                \n\
                usage: ./programName sim_parameter_file\n";

int main(int argc, char *argv[])
{
    if (argc != 2) {
        printf("%s\n", usage);
        exit(1);
    }

    // Parse parameter file
    int i, j, k;
    int nUsers = 0;
    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        printf("File %s not exist\n", argv[1]);
        printf("%s\n", usage);
        exit(1);
    }
    int BUFSIZE = 32768;
    char *buf = calloc(BUFSIZE, sizeof(char));
    // 1, Read number of UEs
    fgets(buf, BUFSIZE, fp);
    sscanf(buf, "%d", &nUsers);
    printf("number of UEs: %d\n", nUsers);

    // 2, Read BS to UE packet loss probabilities
    fgets(buf, BUFSIZE, fp);
    double ploss_B2U[nUsers];
    /*
    double *ploss_B2U = calloc(nUsers, sizeof(double));
    if (ploss_B2U==NULL) {
        fprintf(stderr, "alloc of ploss_B2U failed\n");
        exit(1);
    }
    */
    int SUBSIZE = 20;
    char *subbuf = calloc(SUBSIZE, sizeof(char));
    int offset = 0;
    for (i=0; i<nUsers; i++) {
        memset(subbuf, 0, sizeof(char)*SUBSIZE);
        sscanf(buf, "%s\t%n", subbuf, &offset);
        buf += offset;
        ploss_B2U[i] = (double) atof(subbuf); 
        printf("ploss_B2U[%d]: %.3f\t", i, ploss_B2U[i]);
    }
    printf("\n");

    // 3, Read UE to UE packet loss probabilities
    double ploss_U2U[nUsers*nUsers];
    /*
    double *ploss_U2U = calloc(nUsers*nUsers, sizeof(double));
    if (ploss_U2U==NULL) {
        fprintf(stderr, "alloc of ploss_B2U failed\n");
        exit(1);
    }
    */
    for (i=0; i<nUsers; i++) {
        memset(buf, 0, sizeof(char)*BUFSIZE);
        fgets(buf, BUFSIZE, fp);
        offset = 0;
        for (j=0; j<nUsers; j++) {
            memset(subbuf, 0, sizeof(char)*SUBSIZE);
            sscanf(buf, "%s\t%n", subbuf, &offset);
            buf += offset;
            ploss_U2U[i*nUsers+j] = (double) atof(subbuf); 
            printf("ploss_U2U[%d][%d]: %.3f\t", i, j, ploss_U2U[i*nUsers+j]);
        }
        printf("\n");
    }
    //free(subbuf);
    //free(buf);
    //fclose(fp);

    // Set up simulation
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

    double alpha = 0.05;        // DO NOT change this. The value is in accordance with that in the paper.

    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    char *databuf = malloc(sp.datasize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, databuf, sp.datasize);
    close(rnd);

    struct snc_context *sc;
    /* Create GNC encoding context */
    if ((sc = snc_create_enc_context(databuf, &sp)) == NULL) {
        fprintf(stderr, "Cannot create snc_context.\n");
        exit(1);
    }

    // Create decoders and buffer for each user
    struct snc_decoder **decoders =  malloc(sizeof(struct snc_decoder*) * nUsers);
    struct snc_buffer **buffers = malloc(sizeof(struct snc_buffer*) * nUsers);
    for (i=0; i<nUsers; i++) {
        decoders[i] = snc_create_decoder(&sp, decoder_t);
        buffers[i] = snc_create_buffer(&sp, 1000);  // create a large enough buffer
    }

    clock_t start, stop, dtime = 0;

    struct snc_packet *pkt;    // pointer of coded packet
    int *nuse = calloc(sizeof(int), nUsers);  // count network uses of users
    double *Ecoop = calloc(sizeof(double), nUsers);  // record cooperation energy consumption of each user

    int count1 = 0;

    // Phase I, broadcast the M source packets uncoded to the users. 
    // Using the systematic RLNC would just fulfill this requirement.
    for (k=0; k<M; k++) {
        count1++;
        pkt = snc_generate_packet(sc);
        // BS broadcast to n users
        for (i=0; i<nUsers; i++) {
            nuse[i] += 1;
            if (rand() % 1000 >= ploss_B2U[i] * 1000) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                struct snc_packet *buf_pkt = snc_duplicate_packet(pkt, &sp);
                snc_buffer_packet(buffers[i], buf_pkt);                         // duplicate a copy to save in the buffer, for the purpose of recoding (D2D cooperation)
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
    }
    // check decoder status
    int decoded = 0;
    for (i=0; i<nUsers; i++) {
        if (snc_decoder_finished(decoders[i]))
            decoded += 1;
    }
    int count2 = 0;
    // Phase II: 
    while (decoded != nUsers) {
        count2 += 1;
        // Phase II.a: BS continues to broadcasting to users using orthogonal spectrum
        pkt = snc_generate_packet(sc);
        for (i=0; i<nUsers; i++) {
            if (snc_decoder_finished(decoders[i])) 
                continue;                           // skip completed users

            nuse[i] += 1;
            if (rand() % 1000 >= ploss_B2U[i] * 1000) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                struct snc_packet *buf_pkt = snc_duplicate_packet(pkt, &sp);
                snc_buffer_packet(buffers[i], buf_pkt);              // save a copy in the buffer, for the purpose of D2D comm.
                snc_process_packet(decoders[i], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
        
        // Phase II.b: D2D communication
        // According to Keshthkarjahromi-TNET-2018 Algorithm 2, the user with the most DoF to other users is selected as the 
        // transmitter UE. Here the "most" also needs to take the link erasure probabilities into account.
        int candidate = -1;
        double maxpotential = 0;
        for (i=0; i<nUsers; i++) {
            int u_dof = snc_get_decoder_dof(decoders[i]);
            // calcualte each user's "potential"
            double potential = 0;
            for (j=0; j<nUsers; j++) {
                if (j==i)
                    continue;
                int o_dof = snc_get_decoder_dof(decoders[j]);
                if (o_dof <= u_dof) {
                    potential += 1 - ploss_U2U[nUsers*i+j];         // the term after minus is the loss probability from user i to j
                }
            }
            if (potential > maxpotential) {
                candidate = i;
                maxpotential = potential;
            }
        }
        // User of index 'candidate' is the transmitter user to broadcast to other users via D2D links
        printf("User %d who has %d DoFs is chosen to broadcast to others via D2D.\n", candidate, snc_get_decoder_dof(decoders[candidate]));
        pkt = snc_recode_packet(buffers[candidate], sched_t);
        for (j=0; j<nUsers; j++) {
            if (j == candidate || snc_decoder_finished(decoders[j])) {
                continue;
            }
            if (rand() % 10000 >= ploss_U2U[candidate*nUsers+j] * 10000) {
                struct snc_packet *brd_pkt = snc_duplicate_packet(pkt, &sp);
                struct snc_packet *buf_pkt = snc_duplicate_packet(pkt, &sp);
                snc_buffer_packet(buffers[j], buf_pkt);              // save a copy in the buffer
                snc_process_packet(decoders[j], brd_pkt);
                snc_free_packet(brd_pkt);
            }
        }
        snc_free_packet(pkt);
        
        // Check decoders' status
        decoded = 0;
        for (i=0; i<nUsers; i++) {
            if (snc_decoder_finished(decoders[i]))
                decoded += 1;
        }
    }

    for (i=0; i<nUsers; i++) {
        struct snc_context *dsc = snc_get_enc_context(decoders[i]);
        unsigned char *rec_buf = snc_recover_data(dsc);
        if (memcmp(databuf, rec_buf, sp.datasize) != 0)
            fprintf(stderr, "recovered is NOT identical to original.\n");
        printf("user: %d overhead: %.6f ops: %.6f network-uses: %d\n", i, snc_decode_overhead(decoders[i]), snc_decode_cost(decoders[i]), nuse[i]);
    }
    printf("snum: %d pktsize: %d GF_POWER: %d nBsBroadcast: %d nD2dBroadcast: %d \n", sp.datasize/sp.size_p, sp.size_p, sp.gfpower, count1+count2, count2);
    snc_free_enc_context(sc);
    for (i=0; i<nUsers; i++) 
        snc_free_decoder(decoders[i]);
    return 0;
}

