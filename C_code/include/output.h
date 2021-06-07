#ifndef OUTPUT_HEADER
#define OUTPUT_HEADER

#include "HPO.h"
#include "stack.h"

void fprintf_cluster_PDB(FILE* fp_out,HPO molecules[],stack* cluster,int *mol_no,int* atom_no);
void fprintf_conf_PDB(FILE* fp_out,HPO molecules[],stack clusters[],int number_of_clusters, int threshold, int greater_than_flag);

#endif