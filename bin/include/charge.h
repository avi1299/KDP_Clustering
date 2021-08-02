#ifndef CHARGE_HEADER
#define CHARGE_HEADER

#include "HPO.h"
#include "stack.h"
#include "constants.h"
#include <set>


int fprintf_K_ions_in_cluster(FILE* fp_out, HPO molecules[], K Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag, coordinates boxlength);
void set_mols_of_interest(int mols_of_interest[], stack* cluster);

#endif