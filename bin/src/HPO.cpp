#include "HPO.h"

/*

//Check the minimum distance between 2 points keeping in mind the opposite boundaries of the box are connected
double periodicBoundaryMindist(coordinates point1, coordinates point2, coordinates boxlength)
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

//Return square of the above distace. Saves on computation of sqrt.
double periodicBoundaryMindistSquare(coordinates point1, coordinates point2, coordinates boxlength)
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
*/

//Defining conditions of being connected
//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom
// int connected_molecules(ION *mol1, ION *mol2, coordinates boxlength)
// {
//     return      (strongly_connected_molecules(mol1,mol2,boxlength)||weakly_connected_molecules(mol1,mol2,boxlength));
// }

int connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag)
{
    //If we are comparing the same molecules there is no connection
    if(mol1==mol2)
        return 0;
        //weak connection = 1 | strong connection = 2
    if(strongly_connected_molecules(mol1,mol2,boxlength,PBC_flag))
            return 2;
    if(weakly_connected_molecules(mol1,mol2,boxlength,PBC_flag))
            return 1;
    return 0;
}



int strongly_connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag)
{
    if(PBC_flag)
    return      periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[O2L_1],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[O2L_2],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[O2L_1],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[O2L_2],boxlength)<CUTOFF;
    else
    return      euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[O2L_1])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[O2L_2])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[O2L_1])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[O2L_2])<CUTOFF;   

}

int weakly_connected_molecules(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag)
{
    if(PBC_flag)
    return      periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[OHL_1],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[OHL_2],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[OHL_1],boxlength)<CUTOFF||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[OHL_2],boxlength)<CUTOFF;
    else
    return      euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[OHL_1])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[OHL_2])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[OHL_1])<CUTOFF||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[OHL_2])<CUTOFF;
}

//Here it is at least 1 HOL atom of the first molecule must be in 2.5 nm of one of the second molecule's O2L or OHL atom. Also the OHL atom of the first molecule
//connected to the above HOL is 3.5 nm away from the corresponding O2L or OHL from the above second molecule. 
int connected_molecules_strict(ION *mol1, ION *mol2, coordinates boxlength, int PBC_flag)
{
    if(PBC_flag)
    return      periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[O2L_1],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_1],mol2->posn[O2L_1],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[O2L_2],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_1],mol2->posn[O2L_2],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[OHL_1],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_1],mol2->posn[OHL_1],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_1],mol2->posn[OHL_2],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_1],mol2->posn[OHL_2],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[O2L_1],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_2],mol2->posn[O2L_1],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[O2L_2],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_2],mol2->posn[O2L_2],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[OHL_1],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_2],mol2->posn[OHL_1],boxlength)<CUTOFF_STRICT||
                periodicBoundaryMindistSquare(mol1->posn[HOL_2],mol2->posn[OHL_2],boxlength)<CUTOFF&&periodicBoundaryMindistSquare(mol1->posn[OHL_2],mol2->posn[OHL_2],boxlength)<CUTOFF_STRICT;
    else
    return      euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[O2L_1])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_1],mol2->posn[O2L_1])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[O2L_2])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_1],mol2->posn[O2L_2])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[OHL_1])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_1],mol2->posn[OHL_1])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_1],mol2->posn[OHL_2])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_1],mol2->posn[OHL_2])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[O2L_1])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_2],mol2->posn[O2L_1])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[O2L_2])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_2],mol2->posn[O2L_2])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[OHL_1])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_2],mol2->posn[OHL_1])<CUTOFF_STRICT||
                euclideanDistanceSquare(mol1->posn[HOL_2],mol2->posn[OHL_2])<CUTOFF&&euclideanDistanceSquare(mol1->posn[OHL_2],mol2->posn[OHL_2])<CUTOFF_STRICT;
}

int connected_K_HPO(K *Kmol, ION *HPOmol, coordinates boxlength, int PBC_flag)
{
        if(PBC_flag)
        return  periodicBoundaryMindistSquare(Kmol->posn,HPOmol->posn[O2L_1],boxlength)<CUTOFF_K_O2L||
                periodicBoundaryMindistSquare(Kmol->posn,HPOmol->posn[O2L_2],boxlength)<CUTOFF_K_O2L;
        else
        return  euclideanDistanceSquare(Kmol->posn,HPOmol->posn[O2L_1])<CUTOFF_K_O2L||
                euclideanDistanceSquare(Kmol->posn,HPOmol->posn[O2L_2])<CUTOFF_K_O2L;

}

int connected_SOL_ION(SOL *SOLmol, ION *HPOmol, coordinates boxlength, int PBC_flag)
{
        if(PBC_flag)
        return  periodicBoundaryMindistSquare(SOLmol->posn[HW1],HPOmol->posn[O2L_1], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW1],HPOmol->posn[O2L_2], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW1],HPOmol->posn[OHL_1], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW1],HPOmol->posn[OHL_2], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW2],HPOmol->posn[O2L_1], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW2],HPOmol->posn[O2L_2], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW2],HPOmol->posn[OHL_1], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[HW2],HPOmol->posn[OHL_2], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[OW],HPOmol->posn[HOL_1], boxlength)||
                periodicBoundaryMindistSquare(SOLmol->posn[OW],HPOmol->posn[HOL_2], boxlength); 
        else
        return  euclideanDistanceSquare(SOLmol->posn[HW1],HPOmol->posn[O2L_1])||
                euclideanDistanceSquare(SOLmol->posn[HW1],HPOmol->posn[O2L_2])||
                euclideanDistanceSquare(SOLmol->posn[HW1],HPOmol->posn[OHL_1])||
                euclideanDistanceSquare(SOLmol->posn[HW1],HPOmol->posn[OHL_2])||
                euclideanDistanceSquare(SOLmol->posn[HW2],HPOmol->posn[O2L_1])||
                euclideanDistanceSquare(SOLmol->posn[HW2],HPOmol->posn[O2L_2])||
                euclideanDistanceSquare(SOLmol->posn[HW2],HPOmol->posn[OHL_1])||
                euclideanDistanceSquare(SOLmol->posn[HW2],HPOmol->posn[OHL_2])||
                euclideanDistanceSquare(SOLmol->posn[OW],HPOmol->posn[HOL_1])||
                euclideanDistanceSquare(SOLmol->posn[OW],HPOmol->posn[HOL_2]);                
}


//Funtion to print the details of the ION molecule
void print_ION(ION *mol)
{
        printf("HPO Molecule:\nP  : [%lf,%lf,%lf]\nOHL: [%lf,%lf,%lf]\nHOL: [%lf,%lf,%lf]\nOHL: [%lf,%lf,%lf]\nHOL: [%lf,%lf,%lf]\nO2L: [%lf,%lf,%lf]\nO2L: [%lf,%lf,%lf]\n",
                mol->posn[PL][0],mol->posn[PL][1],mol->posn[PL][2],
                mol->posn[OHL_1][0],mol->posn[OHL_1][1],mol->posn[OHL_1][2],
                mol->posn[HOL_1][0],mol->posn[HOL_1][1],mol->posn[HOL_1][2],
                mol->posn[OHL_2][0],mol->posn[OHL_2][1],mol->posn[OHL_2][2],
                mol->posn[HOL_2][0],mol->posn[HOL_2][1],mol->posn[HOL_2][2],
                mol->posn[O2L_1][0],mol->posn[O2L_1][1],mol->posn[O2L_1][2],
                mol->posn[O2L_2][0],mol->posn[O2L_2][1],mol->posn[O2L_2][2]);
}

//Checks and reports whether the strict definition of connectedness results in the same or lesser connections
void strict_vs_relaxed(ION molecules[],coordinates boxlength,int no_of_molecules, int PBC_flag)
{
    int i,j;
    int strict=0,not_strict=0;
    for(i=0;i<no_of_molecules-1;i++)
        for(j=i+1;j<no_of_molecules;j++)
        {
            //printf("Are molecules %d and %d connected: %d\n",i,j,connected_molecules(molecules[i],molecules[j],boxlength));
            strict+=connected_molecules_strict(&molecules[i],&molecules[j],boxlength, PBC_flag);
            not_strict+=connected_molecules(&molecules[i],&molecules[j],boxlength,PBC_flag);
        }
    printf("Strict : %d  Not Strict: %d\n",strict,not_strict);
}
