#include "graph.h"

//Constructs adjacency list from the array of molecules by checking connectednes between molecules
void adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
    }
    for(i=0;i<no_of_molecules-1;i++)
    {
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(connected_molecules(&molecules[i],&molecules[j],boxlength))
            {
                add_node_given_value(&adjacency_list[i],j);
                add_node_given_value(&adjacency_list[j],i);
            }
        }
    }
}

void adjacency_list_constructor_verbose(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
    }
    for(i=0;i<no_of_molecules-1;i++)
    {
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(connected_molecules(&molecules[i],&molecules[j],boxlength))
            {
                add_node_given_value(&adjacency_list[i],j);
                add_node_given_value(&adjacency_list[j],i);

                //Shows which molecules are connected
                printf("Nodes %d and %d are connected\n",i,j);
                //print_stack(&adjacency_list[i]);
                //print_stack(&adjacency_list[j]);
            }
        }
    }
}