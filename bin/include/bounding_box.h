#ifndef BOUNDING_BOX_HEADER
#define BOUNDING_BOX_HEADER


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define CUTOFF 6.25//Cutoff is 2.5 nm but we are squaring it to save on computation


#define SQR(x) (x)*(x)
//coordinates is a typedef that represents the X,Y,Z values of a point in space
typedef double coordinates[3];


//Check the minimum distance between 2 points keeping in mind the opposite boundaries of the box are connected
double mindist(coordinates point1, coordinates point2, coordinates boxlength);

//Return square of the above distace. Saves on computation of sqrt.
double mindist_square(coordinates point1, coordinates point2, coordinates boxlength);

//Return whether the two atoms are within the cutoff distance i.e. 2.5nm
int connected_or_not(coordinates point1, coordinates point2, coordinates boxlength);

int connected_or_not_NPBC(coordinates point1, coordinates point2);
#endif