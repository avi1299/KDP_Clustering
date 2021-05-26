//#include "HPO.h"
//#include "stack.h"
#include "DFS_clustering.h"
#include "input.h"

//#define LLEN 300
#define NAME 100
#define MAX_M 5000


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

                //Shows which molecules are connected
                printf("Nodes %d and %d are connected\n",i,j);
                //print_stack(&adjacency_list[i]);
                //print_stack(&adjacency_list[j]);
            }
        }
    }
}


int main(int argc,char *argv[])
{
    FILE *fp_in, *fp_out;
    int start_mol_no=-1;
    int no_of_molecules=0;
    int i;
    HPO molecules[MAX_M]; //only for 5000 molecules
    //static cluster_mt cluster[MAX_M];//for all the clusters
    coordinates boxlength;//, coordinate;
    /*------------------------- read the arguments-------------------------*/
    int c;
    while(( c = getopt(argc, argv, "f:o:")) != -1 )
    {
        switch(c)
        {
            case 'f':
                if(( fp_in  =fopen(optarg,"r"))==NULL)
                {
                    printf("cannot open Infile \n"); 
                    exit(0); 
                }
                break;
            case 'o':
                if(( fp_out  =fopen(optarg,"w"))==NULL)
                {
                    printf("cannot open Outfile \n"); 
                    exit(0); 
                }
                break;
        }
    }
    /*------------------------- read the arguments-------------------------*/


    /*---------------------------- read the file --------------------------*/

    //Populate the array after reading PDB
    PDB_reader(fp_in,molecules,boxlength,&no_of_molecules,&start_mol_no);

    //Prints all the molecules
    // printf("%d %d",start_mol_no,end_mol_no);
    // for(i=start_mol_no;i<end_mol_no;i++)
    //     print_HPO(molecules[i]);

    /*---------------------------- read the file --------------------------*/

    //Check if strictness matters
    strict_vs_relaxed(molecules,boxlength,no_of_molecules);



    stack adjacency_list[no_of_molecules];
   // int visited[no_of_molecules]={-1};
    adjacency_list_constructor(molecules,boxlength,no_of_molecules,adjacency_list);

    /*------------------------connectedness calculations-----------------*/
    
    //Printing the adjacency list and connectedness
    for(i=0;i<no_of_molecules;i++)
    {
        printf("Molecule %d\t| Connected to %d molecules | Stack : ",i,adjacency_list[i].length);
        print_stack(&adjacency_list[i]);
    }

    //Statistics of connectedness
    int connectedness[5];
    for(i=0;i<no_of_molecules;i++)
        connectedness[adjacency_list[i].length]++;
    for(i=0;i<5;i++)
        printf("The number of molecules with %d connections is %d\n",i,connectedness[i]);

    /*------------------------connectedness calculations-----------------*/

    int visited[no_of_molecules];
    for(i=0;i<no_of_molecules;i++)
        visited[i]=-1;

    int number_of_clusters;

    dfs(adjacency_list,visited,no_of_molecules,&number_of_clusters);

    /*If the adjacency lists are empty after the operation, it is an indication that
    all the vertices were processed during DFS*/  
    /* 
    for(i=0;i<no_of_molecules;i++)
    {
        printf("%d :",i);
        print_stack(&adjacency_list[i]);
    }
    */

    /*--------------------------clustering calculations------------------*/

    printf("Number of clusters: %d\n",number_of_clusters);


    //Printing the cluster number
    //for(i=0;i<no_of_molecules;i++)
    //    printf("%d molecule belongs to cluster number %d\n",i+start_mol_no,visited[i]);

    stack cluster[number_of_clusters];
    int cluster_max_size=0;
    int cluster_size[no_of_molecules];
    for(i=0;i<no_of_molecules;i++)
        cluster_size[i]=0;

    for(i=0;i<number_of_clusters;i++)
    {
        cluster[i].length=0;
        cluster[i].top=NULL;
    }
    for(i=0;i<no_of_molecules;i++)
        add_node_given_value(&cluster[visited[i]],i);
    for(i=0;i<number_of_clusters;i++)
    {
        printf("Cluster %d\t| Cluster size : %d | Molecules: ",i,cluster[i].length);
        print_stack(&cluster[i]);

        cluster_size[cluster[i].length]++;
        if(cluster[i].length>cluster_max_size)
            cluster_max_size=cluster[i].length;
    }

    for(i=0;i<=cluster_max_size;i++)
        printf("The number of clusters with %d molecules is %d\n",i,cluster_size[i]);



    /*--------------------------clustering calculations------------------*/


    return 0;
}