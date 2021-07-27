#ifndef RING_HEADER
#define RING_HEADER

#include "stack.h"
#include "constants.h"
#include <vector>
#include <set>
#include <utility>

using namespace std;

void makePIDmatrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], 
    int no_of_molecules, int distance_matrix[MAX_MOLECULES][MAX_MOLECULES], 
    vector<vector<vector<set<pair<int,int>>>>> P,
    vector<vector<vector<set<pair<int,int>>>>> P_dash,int strong_flag);

#endif