#ifndef INPUT_HEADER
#define INPUT_HEADER
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "HPO.h"

#define LLEN 300

void PDB_reader(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no);
//void PDB_reader_alternate(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no);
//void PDB_reader_universal(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no);

#endif
