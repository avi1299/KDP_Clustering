#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "HPO.h"
#include "stack.h"

void adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);
void adjacency_list_constructor_verbose(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);
void strongly_adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);
void strongly_adjacency_list_constructor_verbose(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);
void calculate_cluster_coordination_number(stack cluster[],int coordination_no[],double cluster_coordination_no[],int number_of_clusters);
void adjacency_complete(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag);
double strong_connection_ratio(HPO molecules[],coordinates boxlength,int no_of_molecules);

#endif