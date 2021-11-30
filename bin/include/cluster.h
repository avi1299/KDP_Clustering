#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include "stack.h"
#include "ring.h"
#include <vector>

using namespace std;

typedef struct
{
    vector<int> ION_list;
    vector<int> COUTERION_list;
    vector<int> SOL_list;
    ringElementsArray ringElements;
    int charge;
    int hydration;
    double avg_coordination_number;
}t_cluster;

#endif
