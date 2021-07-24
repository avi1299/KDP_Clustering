#ifndef INPUT_HEADER
#define INPUT_HEADER
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include "HPO.h"
#include "gromacs/fileio/xtcio.h"

#define LLEN 300

void XTC_reader(struct t_fileio* fio,FILE* fp_top,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number,real time_to_start);
void PDB_reader(FILE* fp_in,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number);
void TOP_reader(FILE* fp_top,int *no_of_molecules,char *molecule_to_search,char ** atom_list);
void read_mol_order(FILE *fp_top,char *line,char** atom_list);
void skip_comments(FILE *fp_top,char *line);
void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,char ** atom_name_list, int no_of_molecules);

#endif
