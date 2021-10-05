#include "bounding_box.h"

//Returns the minimum distance between coordinates point1 and point2 across the periodic boundary at the boxlength
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
        return euclideanDistance(point3,point2);
}

//Returns the square of minimum distance between coordinates point1 and point2 across the periodic boundary at the boxlength. Better than periodicBoundaryMindist as it saves on computation of sqrt.
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

//Returns the euclidean distance between coordinates point1 and point2
double euclideanDistance(coordinates point1, coordinates point2)
{
        return sqrt(euclideanDistanceSquare(point1,point2));
}


//Returns the square of the euclidean distance between coordinates point1 and point2. Better than euclideanDistance as it saves on computation of sqrt.
double euclideanDistanceSquare(coordinates point1, coordinates point2)
{
        return (SQR(point2[0]-point1[0])+SQR(point2[1]-point1[1])+SQR(point2[2]-point1[2]));
}


// //Return whether the two atoms are within the cutoff distance i.e. 2.5nm
// int connected_or_not(coordinates point1, coordinates point2, coordinates boxlength)
// {
//         int i, j;
//         coordinates point3;
//         for(i=0; i<3; i++)
//                 if(point1[i]-point2[i]>0.5*boxlength[i])
//                         point3[i] = point1[i]-boxlength[i];
//                 else if(point1[i]-point2[i]<-0.5*boxlength[i])
//                         point3[i] = point1[i]+boxlength[i];
//                 else
//                         point3[i]=point1[i];
//         return (SQR(point2[0]-point3[0])+SQR(point2[1]-point3[1])+SQR(point2[2]-point3[2]))<CUTOFF;
// }

// int connected_or_not_NPBC(coordinates point1, coordinates point2)
// {
    
// }