//#include "HPO.h"
//#include "stack.h"
#include "DFS_clustering.h"
#include "input.h"
#include "output.h"
#include "graph.h"
//#include "cluster.h"
#include "constants.h"
#include "charge.h"
#include <omp.h>
#include <string.h>
#include "gromacs/fileio/xtcio.h"
#include <time.h>



int main(int argc,char *argv[])
{
    //Initializing flags, time variables, file pointers and so on
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);
    double overall_cluster_size=0,overall_percentage_clustered=0,overall_percentage_strong=0,overall_cluster_charge=0;
    int cluster_charge;
    double percent_clustered_molecules,strong_ratio;
    FILE *fp_in = NULL, *fp_out=NULL, *fp_top=NULL;
    real time_to_start=0.0;
    int print_every_x_confs=10;
    struct t_fileio* fio = NULL;
    static int verbose_level=0;
    static int check_strict_flag=0;
    static int greater_than_flag=0;
    static int strong_connections_flag=0;
    static int threshold=-1;
    static int threshold_flag=0;
    static int probability_flag=0;
    static int XTC_in_flag=0;
    const char *ext = ".xtc";
    int xlen = strlen(ext);
    int slen;
    /*------------------------- START: read the arguments-------------------------*/
    int c;
    while(( c = getopt(argc, argv, "f:o:v:t:cgs:hpm:r")) != -1 )
    {
        switch(c)
        {
            case 'h':
                printf("Usage: ./test [OPTIONS].. [ARGUMENTS]..\n");
                printf("  -h\t\t: Prints the help information and exits\n");
                printf("  -f <file>\t: Specify the input file\n");
                printf("  -t <file>\t: Specify the .top file in case of XTC input\n");
                printf("  -o <file>\t: Specify the output file\n");
                printf("  -v <int>\t: Specify the verbose level among {1,2,3,4}\n");
                printf("  -c \t\t: Enables checking if stricter connectedness conditions affects number of connections in each configuration\n");
                printf("  -s <int>\t: Specify size of clusters to be outputted in PDB format\n");
                printf("  -g \t\t: Used with -s to include clusters having size greater than or equal to the argument for -s\n");
                printf("  -p \t\t: Prints the percentage of molecule belonging to a cluster rather than the number of molecules\n");
                printf("  -m <int>\t: Prints statisitics every argument number of configurations. Needs verbose levle to be 0\n");
                printf("  -r \t\t: Uses only strong bonds to perform clustering\n");
                exit(0);
            case 'f':
                slen = strlen(optarg);
                if(strcasecmp(optarg + slen - xlen, ext) == 0)
                {
                    XTC_in_flag=1;
                    fio = open_xtc(optarg, "r");
                }
                else if(( fp_in  =fopen(optarg,"r"))==NULL)
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
            case 't':
                if(( fp_top  =fopen(optarg,"r"))==NULL)
                {
                    printf("cannot open Topfile \n"); 
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
            case 'r':
                strong_connections_flag=1;
                break;
            case 'p':
                probability_flag=1;
                break;
            case 'g':
                greater_than_flag=1;
                break;
            case 's':
                threshold_flag=1;
                threshold=atoi(optarg);
                break;
             case 'm':
                print_every_x_confs=atoi(optarg);
                if(print_every_x_confs<1)
                {
                    printf("Please set argument for -m as in int greater than 0\n");
                    exit(0);
                }
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

    FILE *fp_stats=fopen("cluster_statistics.dat","w");
    FILE *fp_cms=fopen("cluster_max_size.dat","w");
    /*------------------------- END: read the arguments-------------------------*/
    int start_mol_no=-1;
    int no_of_molecules=0;
    int i,j;

    static HPO molecules[MAX_M];
    static K Kmolecules[MAX_M];
    static int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES];
    static int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES];
    coordinates boxlength;//, coordinate;
    int conf_number,conf,number_of_K_molecules_in_cluster;
    HPO* mol_start=NULL;
    K* Kmol_start=NULL;
    /*-------------------------START: read the file --------------------------*/
    if(XTC_in_flag==1)
    {
        if(fp_top==NULL)
        {
            printf("Top file not specified\n");
            exit(0);
        }
        XTC_reader(fio,fp_top,molecules,Kmolecules,boxlength,&no_of_molecules,&start_mol_no,&conf_number,time_to_start);

    }
    else
        PDB_reader(fp_in,molecules,Kmolecules,boxlength,&no_of_molecules,&start_mol_no,&conf_number);

    if(verbose_level>=1)
    {
        printf("Number of configurations: %d\n\n",conf_number);
        printf("Number of molecules: %d\n\n",no_of_molecules);
    }
    if (verbose_level>=3)
        printf("BoxLength = %lf %lf %lf \n\n",boxlength[0],boxlength[1],boxlength[2]);

    /*-------------------------END: read the file --------------------------*/
    //Constructing the graph by the means of an adjacency list
    stack adjacency_list[no_of_molecules];
    int coordination_no[no_of_molecules];
    double average_cluster_coordination_no[no_of_molecules];
    double coordination_no_for_cluster_size[no_of_molecules+1];
    int connectedness[MAX_CONNECTIONS+1];
    int visited[no_of_molecules];
    int number_of_clusters;
    stack cluster[no_of_molecules];
    int cluster_max_size=0;
    int cluster_size[no_of_molecules+1];
    double cluster_coordination_statistic[no_of_molecules+1];

    //Printing the start of the PDB file
    if(fp_out!=NULL)
    {
        fprintf(fp_out, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
                boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);
        // printf("CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
        //         boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);     
    }
    
    for(conf=0;conf<conf_number;conf++)
    {
        if (verbose_level>=1)
            printf("Configuration number : %d\n\n",conf+1);
        mol_start=&molecules[conf*no_of_molecules];
        Kmol_start=&Kmolecules[conf*no_of_molecules];
        
        
        //Check if strictness matters
        if(check_strict_flag)
        {
            strict_vs_relaxed(mol_start,boxlength,no_of_molecules);
            printf("\n");
        }
        

        if (verbose_level>=4)
        {
            //Prints all the molecules
            
            for(i=0;i<no_of_molecules;i++)
                print_HPO(&mol_start[i]);
            printf("\n");
        }


        /*------------------------START: connectedness calculations-----------------*/
        
        /*This part is parallelized */
        //Useful when you dont want to claculate the strong ratio
        //adjacency_complete(mol_start,boxlength,no_of_molecules,adjacency_list,(verbose_level>=3),strong_connections_flag);
        adjacency_matrix_populator(mol_start,boxlength,no_of_molecules,adjacency_matrix);
        Kadjacency_matrix_populator(mol_start,Kmol_start,boxlength,no_of_molecules,Kadjacency_matrix);
        adjacency_list_from_matrix(adjacency_matrix,no_of_molecules,adjacency_list,(verbose_level>=3),strong_connections_flag);
    
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

        for(i=0;i<=MAX_CONNECTIONS;i++)
            connectedness[i]=0;    FILE *fp_stats=fopen("cluster_statistics.dat","w");
        FILE *fp_cms=fopen("cluster_max_size.dat","w");
        for(i=0;i<no_of_molecules;i++)
        {
            connectedness[adjacency_list[i].length]++;
            coordination_no[i]=adjacency_list[i].length;
        }
        if(verbose_level>=2)
        {
            //Statistics of connectedness
            for(i=0;i<=MAX_CONNECTIONS;i++)
                printf("The number of molecules with %d connections is %d\n",i,connectedness[i]);
            printf("\n");
        }


        /*-----------------------END: connectedness calculations-----------------*/

        //The visited array is initialized to check if the molecules is visited during DFS.
        //After running dfs it contains the cluster number of the molecules with its index
        for(i=0;i<no_of_molecules;i++)
            visited[i]=-1;

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




        //Printing the cluster number
        //for(i=0;i<no_of_molecules;i++)
        //    printf("%d molecule belongs to cluster number %d\n",i+start_mol_no,visited[i]);

        //Setting the array to 0;
        for(i=0;i<=no_of_molecules;i++)
            cluster_size[i]=0;
            
        for(i=0;i<no_of_molecules;i++)
        {
            average_cluster_coordination_no[i]=0;
            cluster_coordination_statistic[i]=0;
        }
            

        //Initializing the stacks
        for(i=0;i<number_of_clusters;i++)
        {
            cluster[i].length=0;
            cluster[i].top=NULL;
        }

        //Pushing molecule number on the stack corresponding to the cluster number
        for(i=0;i<no_of_molecules;i++)
            add_node_given_value(&cluster[visited[i]],i);

        //Recording cluster_size and cluster_max_size
        cluster_max_size=0;
        for(i=0;i<number_of_clusters;i++)
        {
            cluster_size[cluster[i].length]++;
            if(cluster[i].length>cluster_max_size)
                cluster_max_size=cluster[i].length;
        }

        

        //Calculation average coordination number of cluster
        calculate_cluster_coordination_number(cluster,coordination_no,average_cluster_coordination_no,number_of_clusters);

        for(i=0;i<number_of_clusters;i++)
        {
            cluster_coordination_statistic[cluster[i].length]+=average_cluster_coordination_no[i];
        }
        for(i=2;i<=cluster_max_size;i++)
        {
            if(cluster_size[i]!=0)
                cluster_coordination_statistic[i]/=cluster_size[i];
            else
                cluster_coordination_statistic[i]=0;
        }

        if(verbose_level>=1)
            printf("Number of clusters: %d  with largest size: %d\n\n",number_of_clusters,cluster_max_size);

        if(verbose_level>=3)
        {
            //printf("%d\n",cluster[1].length);
            for(i=0;i<number_of_clusters;i++)
            {
                printf("Cluster %d\t| Cluster size : %d | Coordin. No. : %lf |Molecules: ",i,cluster[i].length,average_cluster_coordination_no[i]);
                print_stack(&cluster[i]);
            }
            printf("\n");
        }


        if (verbose_level>=2)
        {
            if(probability_flag==0)
                for(i=1;i<=cluster_max_size;i++)
                    printf("The number of clusters with %d molecules is %d with avg connections %lf\n",i,cluster_size[i],cluster_coordination_statistic[i]);
            else
                for(i=1;i<=cluster_max_size;i++)
                    printf("Percentage of molecules belonging to a cluster of size %d is %5.3lf\n",i,(double)cluster_size[i]*i/no_of_molecules*100);            
            printf("\n");
        }

        // if((conf%10==0)&&(verbose_level==0))
        // {
        //     printf("Conf: %d | Max Cluster size: %d\n",conf,cluster_max_size);
        // }

        //Adding details to cluster_statisitcs.dat
        for(i=1;i<=no_of_molecules;i++)
            fprintf(fp_stats,"%d ",cluster_size[i]);
        fprintf(fp_stats,"\n");

        fprintf(fp_cms,"%d\n",cluster_max_size);
        // int k;
        // if(conf==0)
        // for(i=0;i<2;i++)
        // for(j=0;j<MAX_MOLECULES;j++)
        // for(k=0;k<MAX_MOLECULES;k++)
        // fprintf(fp_cms,"%d",adjacency_matrix[i][j][k]);

        

        /*--------------------------END: clustering calculations------------------*/

        /*-----------------------START: Output PDB------------------*/
        if(threshold_flag==0)
            threshold=cluster_max_size;
        if(fp_out!=NULL)
        {
            number_of_K_molecules_in_cluster=fprintf_K_ions_in_cluster(fp_out,Kadjacency_matrix,Kmol_start,cluster,no_of_molecules,number_of_clusters,threshold,greater_than_flag);
            fprintf_conf_PDB(fp_out,mol_start,cluster,number_of_clusters,threshold,greater_than_flag,number_of_K_molecules_in_cluster);
        }

        /*-----------------------END: Output PDB------------------*/

        /*-----------------------START: Print Stats------------------*/
        // if(threshold_flag==0)
        //     threshold=cluster_max_size;
        // if(fp_out!=NULL)
        // {
        //     fprintf_conf_PDB(fp_out,mol_start,cluster,number_of_clusters,threshold,greater_than_flag);
        // }

        /*-----------------------END: Print Stats------------------*/
        for(i=0;i<number_of_clusters;i++)
        {
            empty_stack(&cluster[i]);
        }

        overall_cluster_size+=cluster_max_size;
        percent_clustered_molecules= ((double)(no_of_molecules-cluster_size[1]))/no_of_molecules*100;
        overall_percentage_clustered+=percent_clustered_molecules;
        strong_ratio=strong_connection_ratio(adjacency_matrix,no_of_molecules)*100;
        overall_percentage_strong+=strong_ratio;
        
        if(greater_than_flag)
        {
            cluster_charge=0;
            for(i=threshold;i<=cluster_max_size;i++)
            {
                cluster_charge-=i*cluster_size[i];
            }
        }
        else
            cluster_charge=-threshold*cluster_size[threshold];
            
        cluster_charge+=number_of_K_molecules_in_cluster;
        overall_cluster_charge+=cluster_charge;

        if((verbose_level==0)&&(conf%print_every_x_confs==0))
        {
            //double percent_clustered_molecules= ((double)(no_of_molecules-cluster_size[1]))/no_of_molecules*100;
            printf("Conf: %5d | MaxClusterSize: %5d | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %3d\n",conf,cluster_max_size,percent_clustered_molecules,strong_ratio,cluster_charge);
        }

    }
    overall_cluster_size/=conf_number;
    overall_percentage_clustered/=conf_number;
    overall_percentage_strong/=conf_number;
    overall_cluster_charge/=conf_number;
    printf("\nOverall     | MaxClusterSize: %5.2lf | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %5.2lf\n",overall_cluster_size,overall_percentage_clustered,overall_percentage_strong,overall_cluster_charge);

    if(fp_out!=NULL)
        fclose(fp_out);

    fclose(fp_stats);
    fclose(fp_cms);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
 
    printf("The elapsed time is %f seconds\n", elapsed);
    return 0;
}