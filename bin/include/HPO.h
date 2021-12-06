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
#define CUTOFF_SOL_ION 4
#define MAX_CONNECTIONS 6
#define SQR(x) (x)*(x)

//coordinates is a typedef that represents the X,Y,Z values of a point in space
typedef double coordinates[3];

#define HPO_ATOM_COUNT 7

#define ION_ATOM_COUNT HPO_ATOM_COUNT

#define K_ATOM_COUNT 1

#define COUNTERION_ATOM_COUNT K_ATOM_COUNT

#define H2O_ATOM_COUNT 3

#define SOL_ATOM_COUNT H2O_ATOM_COUNT

 /* ------- Defining the molecules that will become the IONS, COUNTERIONS and SOL---------- */

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

typedef struct
{
        coordinates posn[SOL_ATOM_COUNT];
} H20;

#define OW 0
#define HW1 1
#define HW2 2

/* ----------- Giving the molecules an alias of ION and COUNTERION ----------- */

#define ION HPO

#define COUNTERION K

#define SOL H20

//Structure for cluster
// typedef struct
// {
//         stack HPOcluster;
//         stack KMoleculesAroundCluster;
// } clusterStruct;

//Defining conditions of being connected
/**
 * @brief Specifies the requirement for a connection between one ION and another ION :  
 * Checks if at least 1 HOL atom one molecule is within 2.5 nm of the other molecule's O2L or OHL atom. 
 * Returns 0 for no bond, 1 for weak bond and 2 for strong bond
 * 
 * @param mol1 ION*
 * @param mol2 ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag);



/**
 * @brief Specifies the requirement for a strong connection between one ION and another ION : 
 * Checks if at least 1 HOL atom one molecule is within 2.5 nm of one of the other molecule's O2L atom.
 * if(PBC_flag)
 
 * @param mol1 ION*
 * @param mol2 ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int strongly_connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag);


/**
 * @brief Specifies the requirement for a weak connection between one ION and another ION : 
 * Checks if at least 1 HOL atom one molecule is within 2.5 nm of one of the other molecule's OHL atom.
 * 
 * @param mol1 ION*
 * @param mol2 ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int weakly_connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag);
#define PL 0
#define OHL_1 1
#define HOL_1 2

/**
 * @brief Checks if an ION is connected to another ION using a stricter condition :  
 * Checks if at least 1 HOL atom one molecule is within 2.5 nm of the other molecule's O2L or OHL atom. 
 * Also performs a stricter check if the OHL connected to the HOL of the first molecule is within 3.5 nm of the corresponding O2L or OHL from the second molecule.
 * 
 * @param mol1 ION*
 * @param mol2 ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_molecules_strict(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag);

/**
 * @brief Checks if the COUNTERION and ION satisfy the condition for connection:
 * If the K atom is within 3.2 nm of one of the HPO molecule's O2L atom
 * 
 * @param Kmol COUNTERION*
 * @param HPOmol ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_K_HPO(K *Kmol, ION *HPOmol, coordinates boxlength, int PBC_flag);

/**
 * @brief Checks if the SOL and ION satisfy the condition for connection:
 * If either of the H atoms of the H2O molecules is within 2.5 nm of one of the HPO molecule's O2L or OHL atom
 * 
 * @param SOLmol SOL*
 * @param HPOmol ION*
 * @param boxlength coordinates
 * @param PBC_flag int
 * @return int 
 */
int connected_SOL_ION(SOL *SOLmol, ION *HPOmol, coordinates boxlength, int PBC_flag);


/**
 * @brief Print the positions of the atoms of the ION molecule
 * 
 * @param mol ION*
 */
void print_ION(ION *mol);

/**
 * @brief Checks and prints whether the strict definition of connectedness results in the same number of connections
 * 
 * @param molecules ION*
 * @param boxlength coordinates
 * @param no_of_molecules int
 * @param PBC_flag int
 */
void strict_vs_relaxed(ION molecules[],coordinates boxlength,int no_of_molecules, int PBC_flag);

#endif
