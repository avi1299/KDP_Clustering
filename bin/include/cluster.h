#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "stack.h"
#include "ring.h"
#include <vector>

using namespace std;

typedef struct
{
    vector<int> ION_list;
    vector<int> COUTERION_list;
    vector<int> SOL_list;
    int charge;
    double avg_coordination_number;
    ringElementsArray ringElements;
}t_cluster;

#endif
