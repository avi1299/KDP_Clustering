#include "graph.h"
#include "cluster.h"

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

void adjacency_complete(ION molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[],int verbose_flag, int strong_flag, int PCB_flag){
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
                if(strongly_connected_molecules(&molecules[i],&molecules[j],boxlength, PCB_flag)||strongly_connected_molecules(&molecules[j],&molecules[i],boxlength, PCB_flag))
                {
                    add_node_given_value(&adjacency_list[i],j);
                    add_node_given_value(&adjacency_list[j],i);
                    if(verbose_flag)
                        printf("Nodes %5.2lf and %5.2lf are connected\n",i,j);
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
                if(connected_molecules(&molecules[i],&molecules[j],boxlength, PCB_flag)||connected_molecules(&molecules[j],&molecules[i],boxlength, PCB_flag))
                {
                    add_node_given_value(&adjacency_list[i],j);
                    add_node_given_value(&adjacency_list[j],i);
                    if(verbose_flag)
                        printf("Nodes %5.2lf and %5.2lf are connected\n",i,j);
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
                    printf("Nodes %5.2lf and %5.2lf are strongly connected\n",i,j);
            }
        }
    }

}

// void adjacency_list_from_matrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES],int no_of_molecules,vector<t_cluster> cluster,int verbose_flag, int strong_flag){
//     //We are using the undirected graph to create the the adjacency lists so that the direction of the connection is ignored
//     //We are only focusssing on the existence of the connection
//     //int i,j;
//     cluster= new cluster
//     #pragma omp parallel for
//     for(int i=0;i<no_of_molecules;i++)
//     {
//         adjacency_list[i].top=NULL;
//         adjacency_list[i].length=0;        
//     }
//     int strong_level=WEAK;
//     if(strong_flag)
//         strong_level=STRONG;

//     //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
//     #pragma omp parallel for collapse(2) //if (parallelism_enabled)
//     for(int i=0;i<no_of_molecules;i++)
//     {
//         for(int j=0;j<no_of_molecules;j++)
//         {
//             if(adjacency_matrix[UNDIRECTED_GRAPH][i][j]>=strong_level)
//             {
//                 add_node_given_value(&adjacency_list[i],j);
//                 if(verbose_flag)
//                     printf("Nodes %5.2lf and %5.2lf are strongly connected\n",i,j);
//             }
//         }
//     }

// }

//Constructs adjacency matrix from the array of molecules by checking connectednes between molecules
void adjacency_matrix_populator(ION molecules[],coordinates boxlength,int no_of_molecules, int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], int PBC_flag){
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
            adjacency_matrix[DIRECTED_GRAPH][i][j]=connected_molecules(&molecules[i],&molecules[j],boxlength,PBC_flag); 
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

void counterion_adjacency_matrix_populator(ION molecules[], COUNTERION Kmolecules[], coordinates boxlength, int no_of_molecules, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int PBC_flag){
    #pragma omp parallel for collapse(2) //if (parallelism_enabled)
    //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(int i=0;i<no_of_molecules;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            Kadjacency_matrix[i][j]=connected_K_HPO(&Kmolecules[i],&molecules[j],boxlength,PBC_flag); 
        }
    }
}

void SOL_adjacency_matrix_populator(ION molecules[], SOL SOLmolecules[], coordinates boxlength, int no_of_molecules, int no_of_SOL, int SOLadjacency_matrix[MAX_SOL][MAX_MOLECULES], int PBC_flag)
{
    #pragma omp parallel for collapse(2)
    for(int i=0;i<no_of_SOL;i++)
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            SOLadjacency_matrix[i][j]=connected_SOL_ION(&SOLmolecules[i],&molecules[j],boxlength,PBC_flag); 
        }
    }

    // int count=0;
    // for(int j=0;j<no_of_molecules;j++)
    // {
    //     for(int i=0;i<no_of_SOL;i++)
    //     {
    //         //count+=SOLadjacency_matrix[i][j];
    //         if(SOLadjacency_matrix[i][j])
    //             printf("MOL %3d, SOL %5d, Connections: %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf %5.2lf\n",j,i,
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW1],molecules[j].posn[O2L_1], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW1],molecules[j].posn[O2L_2], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW1],molecules[j].posn[OHL_1], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW1],molecules[j].posn[OHL_2], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW2],molecules[j].posn[O2L_1], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW2],molecules[j].posn[O2L_2], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW2],molecules[j].posn[OHL_1], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[HW2],molecules[j].posn[OHL_2], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[OW],molecules[j].posn[HOL_1], boxlength),
    //             periodicBoundaryMindistSquare(SOLmolecules[i].posn[OW],molecules[j].posn[HOL_2], boxlength));
    //             // printf("MOL: %3d, SOL: %5d, Dist: %5.2f, CUTOFF: %5.2lf\n",j,i,periodicBoundaryMindistSquare(SOLmolecules[i].posn[OW],molecules[j].posn[PL], boxlength), CUTOFF);
    //     }
    // }
    // // count/=no_of_molecules;
    // // printf("AVG SOL/MOL = %5.2lf\n",count);

}