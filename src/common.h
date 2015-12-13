/*--------------------- common.h ------------------------
 *  Internal header file 
 *------------------------------------------------------*/
#ifndef COMMON_H
#define COMMON_H
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include "sncEncoder.h"

#ifndef TRACE_LEVEL 
#define TRACE_LEVEL	4
#endif
/* log levels */
#define LOG_ERROR		1
#define LOG_WARNING		2
#define LOG_INFO		3
#define LOG_DEBUG		4
#define LOG_TRACE		5

#define ALIGN(a, b) ((a) % (b) == 0 ? (a)/(b) : (a)/(b) + 1)
#define RESIDUAL(a, b) ((b) * ALIGN((a), (b)) - (a))

#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif
typedef struct node      NBR_node;
typedef struct node_list NBR_nodes;
// node of singly linked list
struct node {
    int data;								
    struct node *next;					
};

struct node_list {
    struct node *first;
    struct node *last;
};

// Bipartitle graph for LDPC code
typedef struct bipartite_graph {
    int		  nleft;
    int		  nright;
    NBR_nodes **l_nbrs_of_r;							// left side neighbours of right
    NBR_nodes **r_nbrs_of_l;							// right side neighbours of left					
} BP_graph;

/**
 * Source packets are grouped into subsets, referred to as 
 * subgenerations in this library to emphasize the 'sub-' concept
 * in sparse code. (Subsets are also referred to as batches, classes, 
 * chunks, segments, or generations in many other network coding literature).
 **/
struct subgeneration {
    int gid;
    int *pktid;							// SIZE_G source packet IDs
};

/**
 * Definition of snc_context
 **/
struct snc_context {
    struct  snc_metainfo   meta;
    struct  subgeneration   **gene; 	// array of pointers each points to a subgeneration.	
    struct  bipartite_graph *graph;
    GF_ELEMENT              **pp;		// Pointers to precoded source packets
};

/* common.c */
int has_item(int array[], int item, int length);
void append_to_list(struct node_list *list, struct node *nd);
int remove_from_list(struct node_list *list, int data);
int exist_in_list(struct node_list *list, int data);
void clear_list(struct node_list *list);
void free_list(struct node_list *list);

/* bipartite.c */
int number_of_checks(int snum, double r); 
int create_bipartite_graph(BP_graph *graph, int nleft, int nright);
void free_bipartite_graph(BP_graph *graph);
#endif /* COMMON_H */
