#include "DFS_clustering.h"

void dfs_util(stack* to_search,stack adjacency_list[],int visited[],int no_of_molecules,int cluster_number)
{
    int c;
    while(to_search->top!=NULL)
    {
        c=pop_and_return_value(to_search);
        visited[c]=cluster_number;
        pop_after_checking_visited(&adjacency_list[c],to_search,visited);
    }
}

void dfs(stack adjacency_list[],int visited[],int no_of_molecules,int *no_of_clusters)
{
    stack to_search;
    to_search.top=NULL;
    int i=0;
    int cluster_number=0;
    while(i<no_of_molecules)
    {
        if(visited[i]==-1)
        {
            visited[i]=cluster_number;
            if(adjacency_list[i].top!=NULL)
            {
                //print_stack(&adjacency_list[i]);
                //add_node_given_value(&to_search,i);
                //print_stack(&to_search);
                pop_after_checking_visited(&adjacency_list[i],&to_search,visited);
                //printf("Popped after checking visited\n");
                //print_stack(&adjacency_list[i]);
                //print_stack(&to_search);
                dfs_util(&to_search, adjacency_list,visited,no_of_molecules,cluster_number);
            }
            cluster_number++;
        }
        i++;
    }
    *no_of_clusters=cluster_number;
}