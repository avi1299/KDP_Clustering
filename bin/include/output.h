#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include "HPO.h"
#include "stack.h"

/**
 * @brief Prints the details of the clusters to the PDB file
 * 
 * @param fp_out FILE*
 * @param molecules ION*
 * @param cluster stack*
 * @param mol_no int*
 * @param atom_no int*
 */
void fprintf_cluster_PDB(FILE* fp_out,ION molecules[],stack* cluster,int *mol_no,int* atom_no);

/** 
 * @brief Prints the details of the configuration to the PDB file
 * 
 * @param fp_out FILE*
 * @param molecules ION*
 * @param clusters stack*
 * @param number_of_clusters int 
 * @param threshold int
 * @param greater_than_flag int 
 * @param HPO_start_no int
 */
void fprintf_conf_PDB(FILE* fp_out,ION molecules[],stack clusters[],int number_of_clusters, int threshold, int greater_than_flag, int HPO_start_no);

#endif