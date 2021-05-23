#include "stack.h"

void dfs(stack adjacency_list[],int visited[],int no_of_molecules);
void dfs_util(stack* to_search,stack adjacency_list[],int visited[],int no_of_molecules,int cluster_number);

