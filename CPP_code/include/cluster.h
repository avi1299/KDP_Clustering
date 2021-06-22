#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "stack.h"

typedef struct
{
    stack list_of_molecules;
    double RoG;
    double avg_coordination_number;
}cluster_wrapper;

#endif
