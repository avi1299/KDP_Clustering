#ifndef BOUNDING_BOX_HEADER
#define BOUNDING_BOX_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

//#define CUTOFF 6.25//Cutoff is 2.5 nm but we are squaring it to save on computation

#define SQR(x) (x)*(x)

//coordinates is a typedef that represents the X,Y,Z values of a point in space
typedef double coordinates[3];

/**
 * @brief Returns the minimum distance between coordinates point1 and point2 across the periodic boundary at the boxlength  
 * 
 * @param point1 coordinates
 * @param point2 coordinates
 * @param boxlength coordinates
 * @return double 
 */
double periodicBoundaryMindist(coordinates point1, coordinates point2, coordinates boxlength);

/**
 * @brief Returns the square of minimum distance between coordinates point1 and point2 across the periodic boundary at the boxlength. Better than periodicBoundaryMindist as it saves on computation of sqrt.
 * 
 * @param point1 coordinates
 * @param point2 coordinates
 * @param boxlength coordinates
 * @return double 
 */
double periodicBoundaryMindistSquare(coordinates point1, coordinates point2, coordinates boxlength);

/**
 * @brief Returns the euclidean distance between coordinates point1 and point2
 * 
 * @param point1 coordinates
 * @param point2 coordinates
 * @return double 
 */
double euclideanDistance(coordinates point1, coordinates point2);

/**
 * @brief Returns the square of the euclidean distance between coordinates point1 and point2. Better than euclideanDistance as it saves on computation of sqrt.
 * 
 * @param point1 coordinates
 * @param point2 coordinates
 * @return double 
 */
double euclideanDistanceSquare(coordinates point1, coordinates point2);

#endif