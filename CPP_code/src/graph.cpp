#include "graph.h"

int parallelism_enabled=1;

//Constructs adjacency list from the array of molecules by checking connectednes between molecules
void adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
    }
    #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        //#pragma omp parallel for
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
    #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        //#pragma omp parallel for
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

void strongly_adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
    }
    #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        //#pragma omp parallel for
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength))
            {
                add_node_given_value(&adjacency_list[i],j);
                add_node_given_value(&adjacency_list[j],i);
            }
        }
    }
}

void strongly_adjacency_list_constructor_verbose(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
    }
    #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        //#pragma omp parallel for
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength))
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

double strong_connection_ratio(HPO molecules[],coordinates boxlength,int no_of_molecules)
{
    int i,j;
    double strong=0,all=0;
    int connection_status=0;
    #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules-1;i++)
    {
        //#pragma omp parallel for
        for(j=i+1;j<no_of_molecules;j++)
        {
            connected_molecules(&molecules[i],&molecules[j],boxlength,&connection_status);
            if(connection_status==2)
            {
                strong++;
                all++;
            }
            else if(connection_status==1)
                all++;
        }
    }
    return strong/all;
}

void adjacency_complete(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag){
    int i,j;
    for(i=0;i<no_of_molecules;i++)
    {
        adjacency_list[i].top=NULL;
        adjacency_list[i].length=0;
        // for (j = 0; i < no_of_molecules; j++){
        //     adjacency_matrix[0][i][j]=0;
        //     adjacency_matrix[1][i][j]=0;
        // }
        
    }
    if(strong_flag)
    {
        #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
        for(i=0;i<no_of_molecules-1;i++)
        {
            //#pragma omp parallel for
            for(j=i+1;j<no_of_molecules;j++)
            {
                if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength))
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
        #pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
        for(i=0;i<no_of_molecules-1;i++)
        {
            //#pragma omp parallel for
            for(j=i+1;j<no_of_molecules;j++)
            {
                if(connected_molecules(&molecules[i],&molecules[j],boxlength))
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