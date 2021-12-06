#ifndef RING_HEADER
#define RING_HEADER

#include "stack.h"
#include "constants.h"
#include <vector>
#include <set>
#include <utility>
#include <limits>
#include <algorithm>

using namespace std;

#define RING_LENGTH_CUTOFF 6

typedef pair<int,int> edge;
typedef vector<edge> path;
typedef vector<path> pathArray;
typedef vector<int> ringElements;
typedef vector<ringElements> ringElementsArray;

typedef struct
{
    int CNum;
    pathArray *P;
    pathArray *P_dash;
} ringCandidate;

struct comp {
    bool operator()(const pair<int,int> &x, const pair<int,int> &y) const {
        if (x.first != y.first) {
            return x.first < y.first;
        }
 
        return x.second < y.second;
    }
};

struct comppathSizeDesc {
    bool operator()(const path &x, const path &y) const { 
        return x.size() > y.size();
    }
};

struct comppathSizeAsc {
    bool operator()(const path &x, const path &y) const { 
        return x.size() < y.size();
    }
};


/**
 * @brief Makes the Path Included Distance Matrix for 
 * 
 * @param adjacency_matrix int *[MAX_MOLECULES][MAX_MOLECULES]
 * @param no_of_molecules int
 * @param D int [MAX_MOLECULES][MAX_MOLECULES]
 * @param P pathArray*[MAX_MOLECULES][MAX_MOLECULES]
 * @param P_dash pathArray*[MAX_MOLECULES][MAX_MOLECULES]
 * @param strong_flag int
 */
void makePIDmatrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag);

/**
 * @brief Searches for the candidate rings from the Path Included Distance matrix
 * 
 * @param CSet vector<ringCandidate> *
 * @param no_of_molecules int
 * @param D int [MAX_MOLECULES][MAX_MOLECULES]
 * @param P pathArray*[MAX_MOLECULES][MAX_MOLECULES]
 * @param P_dash pathArray*[MAX_MOLECULES][MAX_MOLECULES]
 */
void ringCandidateSearch(vector<ringCandidate> *CSet, int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES]);

/**
 * @brief Find the Smallest Set of Smallest Rings from candidate rings
 * 
 * @param CSSSR pathArray *
 * @param CSet vector<ringCandidate> *
 */
void findSSSR(pathArray *CSSSR, vector<ringCandidate> *CSet);

/**
 * @brief Checks if the ring already exists in the smallest set of smallest rings and if it doesn't then it adds the rings to the SSSR
 * 
 * @param CSSSR pathArray *
 * @param ring path*
 */
void pathArrayXORandAdd(pathArray *CSSSR, path *ring);

/**
 * @brief Prints the edge
 * 
 * @param intPair edge*
 */
void printEdge(edge *intPair);

/**
 * @brief Prints the path/ring
 * 
 * @param path path* 
 */
void printPath(path *path);

/**
 * @brief Prints the array of paths/rings
 * 
 * @param patharray pathArray*
 */
void printPathArray(pathArray *patharray);

/**
 * @brief Prints the elements of the rings
 * 
 * @param ring ringElements*
 */
void printRingElements(ringElements *ring);

/**
 * @brief Prints the elements of the array of rings
 * 
 * @param ringarray ringElementsArray *
 */
void printRingElementsArray(ringElementsArray *ringarray);


/**
 * @brief Concatenates each pair of paths in arrays arr1 and arr2 and returns the new pathArray
 * 
 * @param arr1 pathArray*
 * @param arr2 pathArray*
 * @return pathArray* 
 */
pathArray *addPathArray(pathArray* arr1, pathArray* arr2);

/**
 * @brief Appends n paths in arr2 to m paths in arr1 in a pairwise fashion to create n*m paths in arr1
 * 
 * @param arr1 pathArray*
 * @param arr2 pathArray*
 */
void appendPathArray(pathArray* arr1, pathArray* arr2);

/**
 * @brief Removes duplicate edge in the path 
 * 
 * @param x path*
 */
void removeDuplicateEdgePairs(path* x);

/**
 * @brief Removes duplicate rings in the array of rings
 * 
 * @param arr pathArray*
 */
void removeDuplicateRings(pathArray *arr);

/**
 * @brief Checks if the small ring is a part of the big ring
 * 
 * @param small path*
 * @param big path*
 * @return true 
 * @return false 
 */
bool ringSubset(path *small, path *big);

/**
 * @brief Checks if the paths have a common edge
 * 
 * @param small path*
 * @param big path*
 * @return true 
 * @return false 
 */
bool pathIntersection(path *small, path *big);

/**
 * @brief Concatenates each pair of paths in arrays arr1 and arr2 (provided that they do not share a common edge) and returns the new pathArray
 * 
 * @param arr1 pathArray*
 * @param arr2 pathArray*
 * @return pathArray* 
 */
pathArray *addPathArrayWithoutCommonElements(pathArray* arr1, pathArray* arr2);

/**
 * @brief Returns the lements present in teh rings
 * 
 * @param ring 
 * @return ringElements* 
 */
ringElements* ringToElements(path* ring);

void ringDriver(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], int* ION_list, 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag, int verbose_level, vector<ringCandidate> *CSet, 
    pathArray* CSSSR,vector<ringElements> *CSSSR_Elements);


#endif