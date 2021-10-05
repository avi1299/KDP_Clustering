#ifndef CHARGE_HEADER
#define CHARGE_HEADER

#include "HPO.h"
#include "stack.h"
#include "constants.h"
#include <set>

/**
 * @brief Prints the details of the K molecules associated with all the clusters to the fp_out file and returns the total number of K ions printed  
 * 
 * @param fp_out FILE* 
 * @param Kadjacency_matrix int [MAX_MOLECULES][MAX_MOLECULES]
 * @param Kmolecules K*
 * @param clusters stack*
 * @param no_of_molecules int
 * @param number_of_clusters int
 * @param threshold int
 * @param greater_than_flag int 
 * @return int 
 */
int fprintf_K_ions_in_cluster(FILE* fp_out, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , K Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag);

/**
 * @brief Set the mols of interest object
 * 
 * @param mols_of_interest int*
 * @param cluster stack*
 */
void set_mols_of_interest(int mols_of_interest[], stack* cluster);

#endif