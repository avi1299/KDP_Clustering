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
    /*------------------------------ read the arguments-----------------------------------------*/
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
    /*------------------------------ read the arguments-----------------------------------------*/

    //Populate the array after reading PDB
    PDB_reader(fp_in,molecules,boxlength,&no_of_molecules,&start_mol_no);

    //Check if strictness matters
    strict_vs_relaxed(molecules,boxlength,no_of_molecules);

    //int end_mol_no=start_mol_no+no_of_molecules;


    //Prints all the molecules
    // printf("%d %d",start_mol_no,end_mol_no);
    // for(i=start_mol_no;i<end_mol_no;i++)
    //     print_HPO(molecules[i]);

    //int connectedness[MAX_M]={0};

    stack adjacency_list[no_of_molecules];
   // int visited[no_of_molecules]={-1};
    adjacency_list_constructor(molecules,boxlength,no_of_molecules,adjacency_list);
    
    //Printing the adjacency list
    for(i=0;i<no_of_molecules;i++)
    {
        printf("%d :",i);
        print_stack(&adjacency_list[i]);
    }


    int visited[no_of_molecules];
    for(i=0;i<no_of_molecules;i++)
        visited[i]=-1;


    dfs(adjacency_list,visited,no_of_molecules);

    //Printing the cluster number
    for(i=0;i<no_of_molecules;i++)
        printf("%d molecule belongs to cluster number %d\n",i+start_mol_no,visited[i]);

    /*If the adjacency lists are empty after the operation, it is an indication that
    all the vertices were processed during DFS*/  
    /* 
    for(i=0;i<no_of_molecules;i++)
    {
        printf("%d :",i);
        print_stack(&adjacency_list[i]);
    }
    */
     

    return 0;
}