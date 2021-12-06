#ifndef CHARGE_HEADER
#define CHARGE_HEADER

#include "HPO.h"
#include "stack.h"
#include "constants.h"
#include "cluster.h"
#include "string.h"
#include <set>
#include <vector>
#include <assert.h>
/**
 * @brief Prints the details of the K molecules associated with all the clusters to the fp_out file and returns the total number of K ions printed  
 * 
 * @param fp_out FILE* 
 * @param Kadjacency_matrix int [MAX_MOLECULES][MAX_MOLECULES]
 * @param Kmolecules COUNTERION*
 * @param clusters stack*
 * @param no_of_molecules int
 * @param number_of_clusters int
 * @param threshold int
 * @param greater_than_flag int 
 * @return int 
 */
int fprintf_K_ions_in_cluster(FILE* fp_out, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , COUNTERION Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag);

/**
 * @brief Set the mols of interest object
 * 
 * @param mols_of_interest int*
 * @param cluster stack*
 */
void set_mols_of_interest(int mols_of_interest[], stack* cluster);

/**
 * @brief Prints to a file the distribution of IONs associated with COUNTERIONs
 * 
 * @param fp_Kstat FILE*
 * @param Kadjacency_matrix int [MAX_MOLECULES][MAX_MOLECULES]
 * @param no_of_molecules int
 */
void count_counterion_affinity(FILE* fp_Kstat, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int no_of_molecules);

/**
 * @brief Prints the SOL ions
 * 
 * @param fp_out 
 * @param SOLadjacency_matrix 
 * @param SOLmolecules 
 * @param clusters 
 * @param no_of_molecules 
 * @param number_of_clusters 
 * @param no_of_SOL 
 * @param current_atom 
 * @param current_mol 
 * @param threshold 
 * @param greater_than_flag 
 * @return int 
 */
int fprintf_SOL(FILE* fp_out, int SOLadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , SOL SOLmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int no_of_SOL, int current_atom, int current_mol,int threshold, int greater_than_flag);

void add_COUNTERION_to_cluster(int CIONadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int cluster_COUNTERION_matrix[MAX_MOLECULES][MAX_MOLECULES], int no_of_ION, int no_of_CION, t_cluster* clusters, int number_of_clusters);

void add_SOL_to_cluster(int SOL_ION_adjacency_matrix[MAX_SOL][MAX_MOLECULES], int cluster_SOL_matrix[MAX_MOLECULES][MAX_SOL], int no_of_ION, int no_of_SOL, t_cluster* clusters, int number_of_clusters);

#endif