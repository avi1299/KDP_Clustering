#include <stdio.h>
//#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#define CUTOFF 6.25//Cutoff is 2.5 nm but we are squaring it to save on computation
#define CUTOFF_STRICT 12.25//Cutoff for strict method is 3.5nm. Squaring it to save on computation
#define LLEN 300
#define NAME 100
#define MAX_M 5000
#define SQR(x) (x)*(x)
typedef double coordinates[3];
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
typedef struct
{
        int yes_or_no;
        int nodeid[1000];
        int node_number;
} cluster_mt;
double mindist(coordinates point1, coordinates point2, coordinates boxlength)
{
        int i, j;
        coordinates point3;
        for(i=0; i<3; i++)
                if(point1[i]-point2[i]>0.5*boxlength[i])
                        point3[i] = point1[i]-boxlength[i];
                else if(point1[i]-point2[i]<-0.5*boxlength[i])
                        point3[i] = point1[i]+boxlength[i];
                else
                        point3[i]=point1[i];
        return sqrt(SQR(point2[0]-point3[0])+SQR(point2[1]-point3[1])+SQR(point2[2]-point3[2]));
}

double mindist_square(coordinates point1, coordinates point2, coordinates boxlength)
{
        int i, j;
        coordinates point3;
        for(i=0; i<3; i++)
                if(point1[i]-point2[i]>0.5*boxlength[i])
                        point3[i] = point1[i]-boxlength[i];
                else if(point1[i]-point2[i]<-0.5*boxlength[i])
                        point3[i] = point1[i]+boxlength[i];
                else
                        point3[i]=point1[i];
        return (SQR(point2[0]-point3[0])+SQR(point2[1]-point3[1])+SQR(point2[2]-point3[2]));
}


//Return whether the two atoms are within the cutoff distance i.e. 2.5nm
int connected_or_not(coordinates point1, coordinates point2, coordinates boxlength)
{
        int i, j;
        coordinates point3;
        for(i=0; i<3; i++)
                if(point1[i]-point2[i]>0.5*boxlength[i])
                        point3[i] = point1[i]-boxlength[i];
                else if(point1[i]-point2[i]<-0.5*boxlength[i])
                        point3[i] = point1[i]+boxlength[i];
                else
                        point3[i]=point1[i];
        return (SQR(point2[0]-point3[0])+SQR(point2[1]-point3[1])+SQR(point2[2]-point3[2]))<CUTOFF;
}
//Defining conditions of being connected
//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom
int connected_molecules(HPO mol1, HPO mol2, coordinates boxlength)
{
    return      mindist_square(mol1.HOL_1,mol2.O2L_1,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_1,mol2.O2L_2,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_2,mol2.O2L_1,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_2,mol2.O2L_2,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_1,mol2.OHL_1,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_1,mol2.OHL_2,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_2,mol2.OHL_1,boxlength)<CUTOFF||
                mindist_square(mol1.HOL_2,mol2.OHL_2,boxlength)<CUTOFF;
}

//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom. Also the OHL atom of the first molecule
//connected to the above HOL is 3.5 nm away from the corresponding O2L or OHL from the above second molecule. 
int connected_molecules_strict(HPO mol1, HPO mol2, coordinates boxlength)
{
    return      mindist_square(mol1.HOL_1,mol2.O2L_1,boxlength)<CUTOFF&&mindist_square(mol1.OHL_1,mol2.O2L_1,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_1,mol2.O2L_2,boxlength)<CUTOFF&&mindist_square(mol1.OHL_1,mol2.O2L_2,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_1,mol2.OHL_1,boxlength)<CUTOFF&&mindist_square(mol1.OHL_1,mol2.OHL_1,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_1,mol2.OHL_2,boxlength)<CUTOFF&&mindist_square(mol1.OHL_1,mol2.OHL_2,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_2,mol2.O2L_1,boxlength)<CUTOFF&&mindist_square(mol1.OHL_2,mol2.O2L_1,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_2,mol2.O2L_2,boxlength)<CUTOFF&&mindist_square(mol1.OHL_2,mol2.O2L_2,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_2,mol2.OHL_1,boxlength)<CUTOFF&&mindist_square(mol1.OHL_2,mol2.OHL_1,boxlength)<CUTOFF_STRICT||
                mindist_square(mol1.HOL_2,mol2.OHL_2,boxlength)<CUTOFF&&mindist_square(mol1.OHL_2,mol2.OHL_2,boxlength)<CUTOFF_STRICT;
}


//Funtion to print the details of the HPO molecule
void print_HPO(HPO mol)
{
        printf("HPO Molecule:\nP  : [%lf,%lf,%lf]\nOHL: [%lf,%lf,%lf]\nHOL: [%lf,%lf,%lf]\nOHL: [%lf,%lf,%lf]\nHOL: [%lf,%lf,%lf]\nO2L: [%lf,%lf,%lf]\nO2L: [%lf,%lf,%lf]\n",
                mol.P[0],mol.P[1],mol.P[2],
                mol.OHL_1[0],mol.OHL_1[1],mol.OHL_1[2],
                mol.HOL_1[0],mol.HOL_1[1],mol.HOL_1[2],
                mol.OHL_2[0],mol.OHL_2[1],mol.OHL_2[2],
                mol.HOL_2[0],mol.HOL_2[1],mol.HOL_2[2],
                mol.O2L_1[0],mol.O2L_1[1],mol.O2L_1[2],
                mol.O2L_2[0],mol.O2L_2[1],mol.O2L_2[2]);
}
