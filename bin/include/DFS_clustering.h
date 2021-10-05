#ifndef DFS_CLUSTERING_HEADER
#define DFS_CLUSTERING_HEADER

#include "stack.h"

/**
 * @brief Performs Depth-First-Search on the graph to identify connected components of the graph.
 * The visited array is updated with a number corresponding to the node number which indicated with connected component that node belongs to.
 * 
 * @param adjacency_list stack*
 * @param visited int*
 * @param no_of_molecules int
 * @param no_of_clusters int*
 */
void dfs(stack adjacency_list[],int visited[],int no_of_molecules, int *no_of_clusters);

/**
 * @brief A utility function for the dfs function.
 * 
 * @param to_search stack*
 * @param adjacency_list stack*
 * @param visited int*
 * @param no_of_molecules int 
 * @param cluster_number int
 */
void dfs_util(stack* to_search,stack adjacency_list[],int visited[],int no_of_molecules,int cluster_number);

#endif

