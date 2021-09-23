#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "stack.h"
#include "ring.h"

typedef struct
{
    stack clusterElements;
    int charge;
    double avg_coordination_number;
    ringElementsArray ringElements;
}cluster_wrap;

#endif
