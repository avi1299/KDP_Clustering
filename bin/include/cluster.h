#ifndef CLUSTER_HEADER
#define CLUSTER_HEADER

#include "stack.h"
#include "ring.h"
#include <vector>

using namespace std;

typedef struct
{
    int* ION_list;
    int ION_list_size;
    int* COUTERION_list;
    int COUNTERION_list_size;
    int* SOL_list;
    int SOL_list_size;
    ringElementsArray ringElements;
    //double avg_coordination_number;
}t_cluster;

#endif
