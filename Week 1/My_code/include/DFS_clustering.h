#include "stack.h"

void dfs(stack adjacency_list[],int visited[],int no_of_molecules, int *no_of_clusters);
void dfs_util(stack* to_search,stack adjacency_list[],int visited[],int no_of_molecules,int cluster_number);

