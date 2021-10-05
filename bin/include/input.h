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

/**
 * @brief Holds the name and quantity of molecules in the TOP file
 * 
 */
typedef struct
{
    string name;
    int quantity;
} orderOfMols;

typedef vector<orderOfMols> orderOfMolsArray;

/**
 * @brief Stores the details of the residues present in the TOP file.
 * 
 */
typedef struct
{
    string moleculeName;
    vector<string> atomName;
    //TODO: Implement Center of Gravity, charge fields in atoms subsection
    //TODO: Implemet bonds, angles and dihedrals
}moleculeInfo;

typedef vector<moleculeInfo> moleculeInfoList; 


/**
 * @brief Reads the XTC file and stores the positions of the molecules
 * 
 * @param fio FILE*
 * @param fp_top FILE*
 * @param molecules HPO*
 * @param Kmolecules K*
 * @param boxlength coordinates
 * @param no_of_molecules int*
 * @param start_mol_no int*
 * @param conf_number int*
 * @param time_to_start real
 */
void XTC_reader(struct t_fileio* fio,FILE* fp_top,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number,real time_to_start);

/**
 * @brief Reads the PDB file and stores the positions of the molecules
 * 
 * @param fp_in FILE*
 * @param molecules HPO*
 * @param Kmolecules K*
 * @param boxlength coordinates
 * @param no_of_molecules int*
 * @param start_mol_no int*
 * @param conf_number int*
 */
void PDB_reader(FILE* fp_in,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number);

/**
 * @brief [Depreceated] Reads the TOP file and stores the number of HPO and K molecules
 * 
 * @param fp_top FILE*
 * @param no_of_molecules int*
 * @param molecule_to_search char*
 * @param atom_list char**
 * @param order orderOfMolsArray*
 */
void TOP_reader(FILE* fp_top,int *no_of_molecules,char *molecule_to_search,char ** atom_list, orderOfMolsArray *order);

/**
 * @brief Parses the TOP file to capture information about the type of residues involved along with their orders and quantities
 * @param fp_top FILE*
 * @param molInfo moleculeInfoList *
 * @param order orderOfMolsArray *
 */
void TOP_parser(FILE* fp_top, moleculeInfoList *molInfo ,orderOfMolsArray *order);

/**
 * @brief [Depreceated] Reads the order of atoms in the molecule
 * 
 * @param fp_top FILE*
 * @param line char*
 * @param atom_list char**
 */
void read_mol_order(FILE *fp_top,char *line,char** atom_list);

/**
 * @brief [Depreceated] Skips commnets in the TOP file
 * 
 * @param fp_top FILE*
 * @param line char*
 */
void skip_comments(FILE *fp_top,char *line);
//void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,char ** atom_name_list, int no_of_molecules);

/**
 * @brief Reads the data from the XTC file into the arrays storeing the positions of the molecules
 * 
 * @param molecules HPO*
 * @param Kmolecules K*
 * @param x rvec*
 * @param atom_index int* 
 * @param no_of_molecules int 
 */
void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,int * atom_index, int no_of_molecules);

/**
 * @brief Converts an array of strings into their corresponding indices
 * 
 * @param atom_name_list char**
 * @param atom_index int*
 */
void atom_name_to_index(char ** atom_name_list, int* atom_index);

/**
 * @brief Converts a vector of strings into their corresponding indices
 * 
 * @param atom_name_list vector<string>
 * @param atom_index vector<int>*
 */
void atom_name_to_index(vector<string> atom_name_list, vector<int> *atom_index);


#endif
