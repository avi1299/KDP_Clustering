#include "bounding_box.h"

//Check the minimum distance between 2 points keeping in mind the opposite boundaries of the box are connected
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

//Return square of the above distace. Saves on computation of sqrt.
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