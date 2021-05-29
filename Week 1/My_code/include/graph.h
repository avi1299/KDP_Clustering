#ifndef GRAPH_HEADER
#define GRAPH_HEADER

#include "HPO.h"
#include "stack.h"

void adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);
void adjacency_list_constructor_verbose(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[]);


#endif