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




void makePIDmatrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag);

void ringCandidateSearch(vector<ringCandidate> *CSet, int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES]);

void findSSSR(pathArray *CSSSR, vector<ringCandidate> *CSet);
void pathArrayXORandAdd(pathArray *CSSSR, path *ring);

void printEdge(edge *intPair);
void printPath(path *path);
void printPathArray(pathArray *patharray);
void printRingElements(ringElements *ring);
void printRingElementsArray(ringElementsArray *ringarray);

pathArray *addPathArray(pathArray* arr1, pathArray* arr2);
void appendPathArray(pathArray* arr1, pathArray* arr2);
void removeDuplicateEdgePairs(path* x);
void removeDuplicateRings(pathArray *arr);
bool ringSubset(path *small, path *big);
bool pathIntersection(path *small, path *big);
pathArray *addPathArrayWithoutCommonElements(pathArray* arr1, pathArray* arr2);
ringElements* ringToElements(path* ring);


#endif