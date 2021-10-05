#ifndef HPO_HEADER
#define HPO_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "bounding_box.h"

#define CUTOFF 6.25//Cutoff is 2.5 nm but we are squaring it to save on computation
#define CUTOFF_STRICT 12.25//Cutoff for strict method is 3.5nm. Squaring it to save on computation
#define CUTOFF_K_O2L 10.24//Cutoff for K-O2L interaction 3.2nm
#define MAX_CONNECTIONS 6
#define SQR(x) (x)*(x)

//coordinates is a typedef that represents the X,Y,Z values of a point in space
typedef double coordinates[3];

#define HPO_ATOM_COUNT 7

typedef struct
{
        coordinates posn[HPO_ATOM_COUNT];
} HPO;

#define PL 0
#define OHL_1 1
#define HOL_1 2
#define OHL_2 3
#define HOL_2 4
#define O2L_1 5
#define O2L_2 6

typedef struct
{
        coordinates posn;
} K;



//Structure for cluster
// typedef struct
// {
//         stack HPOcluster;
//         stack KMoleculesAroundCluster;
// } clusterStruct;

//Defining conditions of being connected
/**
 * @brief Checks if at least 1 HOL atom one molecule is within 2.5 nm of the other molecule's O2L or OHL atom. 
 * Returns 0 for no bond, 1 for weak bond and 2 for strong bond
 * 
 * @param mol1 HPO*
 * @param mol2 HPO*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_molecules(HPO *mol1, HPO *mol2, coordinates boxlength, int PBC_flag);



/**
 * @brief Checks if at least 1 HOL atom one molecule is within 2.5 nm of one of the other molecule's O2L atom.
 * 
 * @param mol1 HPO*
 * @param mol2 HPO*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int strongly_connected_molecules(HPO *mol1, HPO *mol2, coordinates boxlength, int PBC_flag);


/**
 * @brief Checks if at least 1 HOL atom one molecule is within 2.5 nm of one of the other molecule's OHL atom.
 * 
 * @param mol1 HPO*
 * @param mol2 HPO*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int weakly_connected_molecules(HPO *mol1, HPO *mol2, coordinates boxlength, int PBC_flag);


/**
 * @brief Checks if at least 1 HOL atom one molecule is within 2.5 nm of the other molecule's O2L or OHL atom. 
 * Also performs a stricter check if the OHL connected to the HOL of the first molecule is within 3.5 nm of the corresponding O2L or OHL from the second molecule.
 * 
 * @param mol1 HPO*
 * @param mol2 HPO*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_molecules_strict(HPO *mol1, HPO *mol2, coordinates boxlength, int PBC_flag);

/**
 * @brief Checks if the K atom is within 3.2 nm of one of the HPO molecule's O2L atom
 * 
 * @param Kmol K*
 * @param HPOmol HPO*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_K_HPO(K *Kmol, HPO *HPOmol, coordinates boxlength, int PBC_flag);


/**
 * @brief Print the positions of the atoms of the HPO molecule
 * 
 * @param mol HPO*
 */
void print_HPO(HPO *mol);

/**
 * @brief Checks and prints whether the strict definition of connectedness results in the same number of connections
 * 
 * @param molecules HPO*
 * @param boxlength coordinates
 * @param no_of_molecules int
 * @param PBC_flag int
 */
void strict_vs_relaxed(HPO molecules[],coordinates boxlength,int no_of_molecules, int PBC_flag);

#endif
