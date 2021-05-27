#ifndef HPO_HEADER
#define HPO_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define CUTOFF 6.25//Cutoff is 2.5 nm but we are squaring it to save on computation
#define CUTOFF_STRICT 12.25//Cutoff for strict method is 3.5nm. Squaring it to save on computation
#define SQR(x) (x)*(x)

//coordinates is a typedef that represents the X,Y,Z values of a point in space
typedef double coordinates[3];

//Structure of the HPO molecule
typedef struct
{
        coordinates P;
        coordinates OHL_1;
        coordinates HOL_1;
        coordinates OHL_2;
        coordinates HOL_2;
        coordinates O2L_1;
        coordinates O2L_2;
} HPO;

//Structure for cluster
typedef struct
{
        int yes_or_no;
        int nodeid[1000];
        int node_number;
} cluster_mt;

//Check the minimum distance between 2 points keeping in mind the opposite boundaries of the box are connected
double mindist(coordinates point1, coordinates point2, coordinates boxlength);

//Return square of the above distace. Saves on computation of sqrt.
double mindist_square(coordinates point1, coordinates point2, coordinates boxlength);

//Return whether the two atoms are within the cutoff distance i.e. 2.5nm
int connected_or_not(coordinates point1, coordinates point2, coordinates boxlength);

//Defining conditions of being connected
//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom
int connected_molecules(HPO *mol1, HPO *mol2, coordinates boxlength);

//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom. Also the OHL atom of the first molecule
//connected to the above HOL is 3.5 nm away from the corresponding O2L or OHL from the above second molecule. 
int connected_molecules_strict(HPO *mol1, HPO *mol2, coordinates boxlength);


//Funtion to print the details of the HPO molecule
void print_HPO(HPO *mol);

//Checks and reports whether the strict definition of connectedness results in the same or lesser connections
void strict_vs_relaxed(HPO molecules[],coordinates boxlength,int no_of_molecules);

#endif
