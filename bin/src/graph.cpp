#include "graph.h"

//int parallelism_enabled=1;

void calculate_cluster_coordination_number(stack cluster[],int coordination_no[],double cluster_coordination_no[],int number_of_clusters)
{
    int i;
    node* temp;
    double sum;
    for(i=0;i<number_of_clusters;i++)
    {
        sum=0;
        temp=cluster[i].top;
        while(temp!=NULL)
        {
            sum+=coordination_no[temp->data];
            temp=temp->next;
        }
        cluster_coordination_no[i]=sum/cluster[i].length;
    }
}

double strong_connection_ratio(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules)
{
    //We are comparing on the ration of strong connections on the directed graph.
    //Doing it on our undirected graph leads to inflated numbers
    //int i,j;
    int strong[no_of_molecules]={0},all[no_of_molecules]={0};
    //Enabling parallellism seems to be messing with count due to race conditions
    //#pragma omp parallel for schedule(static, 1) private(j)  if (parallelism_enabled)
    #pragma omp parallel for //if (parallelism_enabled)
    for(int i=0;i<no_of_molecules;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            //strong connection
            if(adjacency_matrix[DIRECTED_GRAPH][i][j]==STRONG)
            {
                strong[i]++;
                all[i]++;
            }
            else if(adjacency_matrix[DIRECTED_GRAPH][i][j]==WEAK)//weak connection
                all[i]++;
        }
    }
    for(int i=1;i<no_of_molecules;i++)
    {
        strong[0]+=strong[i];
        all[0]+=all[i];
    }
    //printf("Strong: %5d | All: %5d\n",strong,all);
    return (double)strong[0]/all[0];
}

void adjacency_complete(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag){
    int i,j;
    #pragma omp parallel for
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
        
    }
    if(strong_flag)
    {
        #pragma omp parallel for schedule(static, 1) private(j) //if (parallelism_enabled)
        for(i=0;i<no_of_molecules-1;i++)
        {
            for(j=i+1;j<no_of_molecules;j++)
            {
                //if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength)||)
                if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength)||strongly_connected_molecules(&molecules[j],&molecules[i],boxlength))
                {
                    add_node_given_value(&adjacency_list[i],j);
                    add_node_given_value(&adjacency_list[j],i);
                    if(verbose_flag)
                        printf("Nodes %d and %d are connected\n",i,j);
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for schedule(static, 1) private(j) //if (parallelism_enabled)
        for(i=0;i<no_of_molecules-1;i++)
        {
            for(j=i+1;j<no_of_molecules;j++)
            {
                //if(connected_molecules(&molecules[i],&molecules[j],boxlength))
                if(connected_molecules(&molecules[i],&molecules[j],boxlength)||connected_molecules(&molecules[j],&molecules[i],boxlength))
                {
                    add_node_given_value(&adjacency_list[i],j);
                    add_node_given_value(&adjacency_list[j],i);
                    if(verbose_flag)
                        printf("Nodes %d and %d are connected\n",i,j);
                }
            }
        }
    }
}

//Constructs adjacency matrix from the array of molecules by checking connectednes between molecules
void adjacency_list_from_matrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag){
    //We are using the undirected graph to create the the adjacency lists so that the direction of the connection is ignored
    //We are only focusssing on the existence of the connection
    //int i,j;
    #pragma omp parallel for
    for(int i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;        
    }
    int strong_level=WEAK;
    if(strong_flag)
        strong_level=STRONG;

    //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    #pragma omp parallel for collapse(2) //if (parallelism_enabled)
    for(int i=0;i<no_of_molecules;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            if(adjacency_matrix[UNDIRECTED_GRAPH][i][j]>=strong_level)
            {
                add_node_given_value(&adjacency_list[i],j);
                if(verbose_flag)
                    printf("Nodes %d and %d are strongly connected\n",i,j);
            }
        }
    }

}

//Constructs adjacency matrix from the array of molecules by checking connectednes between molecules
void adjacency_matrix_populator(HPO molecules[],coordinates boxlength,int no_of_molecules, int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES]){
    #pragma omp parallel for collapse(2) //if (parallelism_enabled)
    //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(int i=0;i<no_of_molecules;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            // if(i!=j)
            //     adjacency_matrix[DIRECTED_GRAPH][i][j]=connected_molecules(&molecules[i],&molecules[j],boxlength);    
            // else
            //     adjacency_matrix[DIRECTED_GRAPH][i][j]=0;
            adjacency_matrix[DIRECTED_GRAPH][i][j]=connected_molecules(&molecules[i],&molecules[j],boxlength); 
        }
    }

    //Creating the undirected graph.
    int i,j;
    //#pragma omp parallel for collapse(2) if (parallelism_enabled)
    #pragma omp parallel for schedule(static, 1) private(j) //if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(adjacency_matrix[DIRECTED_GRAPH][i][j]>adjacency_matrix[DIRECTED_GRAPH][j][i]){
                adjacency_matrix[UNDIRECTED_GRAPH][i][j]=adjacency_matrix[DIRECTED_GRAPH][i][j];
                adjacency_matrix[UNDIRECTED_GRAPH][j][i]=adjacency_matrix[DIRECTED_GRAPH][i][j];
            }
            else{
                adjacency_matrix[UNDIRECTED_GRAPH][i][j]=adjacency_matrix[DIRECTED_GRAPH][j][i];
                adjacency_matrix[UNDIRECTED_GRAPH][j][i]=adjacency_matrix[DIRECTED_GRAPH][j][i];
            }
                
        }
    }

}

void Kadjacency_matrix_populator(HPO molecules[], K Kmolecules[], coordinates boxlength, int no_of_molecules, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES]){
    #pragma omp parallel for collapse(2) //if (parallelism_enabled)
    //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(int i=0;i<no_of_molecules;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            Kadjacency_matrix[i][j]=connected_K_HPO(&Kmolecules[i],&molecules[j],boxlength); 
        }
    }
}