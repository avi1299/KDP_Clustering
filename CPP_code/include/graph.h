#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "HPO.h"
#include "stack.h"
#include "constants.h"

void calculate_cluster_coordination_number(stack cluster[],int coordination_no[],double cluster_coordination_no[],int number_of_clusters);
void adjacency_complete(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag);
void adjacency_matrix_populator(HPO molecules[],coordinates boxlength,int no_of_molecules, int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES]);
void adjacency_list_from_matrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag);
double strong_connection_ratio(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules);


#endif