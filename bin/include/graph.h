#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "HPO.h"
#include "stack.h"
#include "constants.h"

/**
 * @brief Calculates the average coordination number(number of connections per molecule) of each cluster 
 * 
 * @param cluster stack*
 * @param coordination_no int*
 * @param cluster_coordination_no double* 
 * @param number_of_clusters int
 */
void calculate_cluster_coordination_number(stack cluster[],int coordination_no[],double cluster_coordination_no[],int number_of_clusters);

/**
 * @brief [Depreceated] Constructs the adjacency list from the graph
 * 
 * @param molecules HPO*
 * @param boxlength coordinates
 * @param no_of_molecules int
 * @param adjacency_list stack*
 * @param verbose_flag int
 * @param strong_flag int
 * @param PCB_flag int
 */
void adjacency_complete(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag, int PCB_flag);

/**
 * @brief Constructs the adjacency matrix from the graph
 * 
 * @param molecules HPO*
 * @param boxlength coordinates
 * @param no_of_molecules int
 * @param adjacency_matrix int* [MAX_MOLECULES][MAX_MOLECULES]
 * @param PBC_flag int
 */
void adjacency_matrix_populator(HPO molecules[],coordinates boxlength,int no_of_molecules, int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int PBC_flag);

/**
 * @brief Constructs the adjacency list from the adjcency matrix
 * 
 * @param adjacency_matrix int* [MAX_MOLECULES][MAX_MOLECULES]
 * @param no_of_molecules int
 * @param adjacency_list stack*
 * @param verbose_flag int
 * @param strong_flag int
 */
void adjacency_list_from_matrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag);

/**
 * @brief Returns the ratio of strong connections to all connections
 * 
 * @param adjacency_matrix int* [MAX_MOLECULES][MAX_MOLECULES]
 * @param no_of_molecules int
 * @return double 
 */
double strong_connection_ratio(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules);

/**
 * @brief Constructs an adjacency matrix which indicates which K molecules are connected to which HPO molecules
 * 
 * @param molecules HPO*
 * @param Kmolecules K*
 * @param boxlength coordinates
 * @param no_of_molecules int
 * @param Kadjacency_matrix int[MAX_MOLECULES][MAX_MOLECULES]
 * @param PBC_flag int
 */
void Kadjacency_matrix_populator(HPO molecules[], K Kmolecules[], coordinates boxlength, int no_of_molecules, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int PBC_flag);


#endif