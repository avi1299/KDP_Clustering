//#include "HPO.h"
//#include "stack.h"
#include "DFS_clustering.h"
#include "input.h"
#include "output.h"

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


int main(int argc,char *argv[])
{
    FILE *fp_in, *fp_out;
    int start_mol_no=-1;
    int no_of_molecules=0;
    int i;
    int verbose_level=0;
    int check_strict_flag=0;
    int greater_than_flag=0;
    int threshold=-1;
    HPO molecules[MAX_M]; //only for 5000 molecules
    //static cluster_mt cluster[MAX_M];//for all the clusters
    coordinates boxlength;//, coordinate;
    /*------------------------- START: read the arguments-------------------------*/
    int c;
    while(( c = getopt(argc, argv, "f:o:v:cgs:")) != -1 )
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
            case 'v':
                verbose_level=atoi(optarg);
                if((verbose_level<1)||(verbose_level>4))
                {
                    printf("Unknown Verbose level = %s\n",optarg);
                    printf("Verbose level 1 -> Number of clusters and max size\n");
                    printf("Verbose level 2 -> Statistics of connectedness and molecules\n");
                    printf("Verbose level 3 -> Complete information of clusters and molecule neighbours\n");
                    printf("Verbose level 4 -> Verbose level 3 and molecular postions of all molecules\n");
                    exit(0);
                }
                break;
            case 'c':
                check_strict_flag=1;
                break;
            case 'g':
                greater_than_flag=1;
                break;
            case 's':
                threshold=atoi(optarg);
                break;
            case '?':
                if (optopt=='v')
                {
                    printf("Verbose level 1 -> Number of clusters and max size\n");
                    printf("Verbose level 2 -> Statistics of connectedness and molecules\n");
                    printf("Verbose level 3 -> Complete information of clusters and molecule neighbours\n");
                    printf("Verbose level 4 -> Verbose level 3 and molecular postions of all molecules\n");
                }
                exit(0);

        }
    }
    /*------------------------- END: read the arguments-------------------------*/


    /*-------------------------START: read the file --------------------------*/

    PDB_reader(fp_in,molecules,boxlength,&no_of_molecules,&start_mol_no);

    //Check if strictness matters
    if(check_strict_flag)
    {
        strict_vs_relaxed(molecules,boxlength,no_of_molecules);
        printf("\n");
    }

    if(verbose_level>=1)
        printf("Number of molecules: %d\n\n",no_of_molecules);
    if (verbose_level>=3)
        printf("BoxLength = %lf %lf %lf \n\n",boxlength[0],boxlength[1],boxlength[2]);
    if (verbose_level>=4)
    {
        //Prints all the molecules
        for(i=0;i<no_of_molecules;i++)
            print_HPO(&molecules[i]);
        printf("\n");
    }

    /*-------------------------END: read the file --------------------------*/

    //Constructing the graph by the means of an adjacency list
    stack adjacency_list[no_of_molecules];
    if(verbose_level>=3)
    {
        adjacency_list_constructor_verbose(molecules,boxlength,no_of_molecules,adjacency_list);
        printf("\n");
    }
    else
        adjacency_list_constructor(molecules,boxlength,no_of_molecules,adjacency_list);

    /*------------------------START: connectedness calculations-----------------*/
    
    //Printing the adjacency list and connectedness
    if(verbose_level>=3)
    {
        for(i=0;i<no_of_molecules;i++)
        {
            printf("Molecule %d\t| Connected to %d molecules | Stack : ",i,adjacency_list[i].length);
            print_stack(&adjacency_list[i]);
        }
        printf("\n");
    }


    //Statistics of connectedness
    int connectedness[5];
    for(i=0;i<no_of_molecules;i++)
        connectedness[adjacency_list[i].length]++;

    if(verbose_level>=2)
    {
        for(i=0;i<5;i++)
            printf("The number of molecules with %d connections is %d\n",i,connectedness[i]);
        printf("\n");
    }


    /*-----------------------END: connectedness calculations-----------------*/

    //The visited array is initialized to check if the molecules is visited during DFS.
    //After running dfs it contains the cluster number of the molecules with its index
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

    /*-----------------------START: clustering calculations------------------*/

    if(verbose_level>=1)
    {
        printf("Number of clusters: %d\n",number_of_clusters);
        printf("\n");
    }


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
        cluster_size[cluster[i].length]++;
        if(cluster[i].length>cluster_max_size)
            cluster_max_size=cluster[i].length;
    }

    if(verbose_level>=3)
    {
        for(i=0;i<number_of_clusters;i++)
        {
            printf("Cluster %d\t| Cluster size : %d | Molecules: ",i,cluster[i].length);
            print_stack(&cluster[i]);
        }
        printf("\n");
    }


    if (verbose_level>=2)
    {
        for(i=1;i<=cluster_max_size;i++)
        printf("The number of clusters with %d molecules is %d\n",i,cluster_size[i]);
        printf("\n");
    }

    /*--------------------------END: clustering calculations------------------*/

    /*-----------------------START: Output PDB------------------*/
    if(threshold==-1)
        threshold=cluster_max_size;
    if(fp_out!=NULL)
    {
        fprintf(fp_out, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
                boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);
        fprintf_conf_PDB(fp_out,molecules,cluster,number_of_clusters,threshold,greater_than_flag);
        fclose(fp_out);
    }

    /*-----------------------END: Output PDB------------------*/


    return 0;
}