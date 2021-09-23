#ifndef INPUT_HEADER
#define INPUT_HEADER
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <vector>

#include "HPO.h"
#include "gromacs/fileio/xtcio.h"

#define LLEN 300
using namespace std;

typedef struct
{
    string name;
    int quantity;
} orderOfMols;

typedef vector<orderOfMols> orderOfMolsArray;

typedef struct
{
    string moleculeName;
    vector<string> atomName;
    //TODO: Implement Center of Gravity, charge fields in atoms subsection
    //TODO: Implemet bonds, angles and dihedrals
}moleculeInfo;

typedef vector<moleculeInfo> moleculeInfoList; 



void XTC_reader(struct t_fileio* fio,FILE* fp_top,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number,real time_to_start);
void PDB_reader(FILE* fp_in,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number);
void TOP_reader(FILE* fp_top,int *no_of_molecules,char *molecule_to_search,char ** atom_list, orderOfMolsArray *order);
void TOP_parser(FILE* fp_top, moleculeInfoList *molInfo ,orderOfMolsArray *order);
void read_mol_order(FILE *fp_top,char *line,char** atom_list);
void skip_comments(FILE *fp_top,char *line);
//void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,char ** atom_name_list, int no_of_molecules);
void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,int * atom_index, int no_of_molecules);
void atom_name_to_index(char ** atom_name_list, int* atom_index);
void atom_name_to_index(vector<string> atom_name_list, vector<int> *atom_index);


#endif
