//#include "HPO.h"
//#include "stack.h"
#include "DFS_clustering.h"
#include "input.h"
#include "output.h"
#include "graph.h"
#include "ring.h"
//#include "cluster.h"
#include "constants.h"
#include "charge.h"
#include <omp.h>
#include <string.h>
#include "gromacs/fileio/xtcio.h"
#include <time.h>

using namespace std;



int main(int argc,char *argv[])
{
    //printf("start\n");
    //Initializing flags, time variables, file pointers and so on
    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);

    double overall_cluster_size=0;
    double overall_percentage_clustered=0;
    double overall_percentage_strong=0;
    double overall_cluster_charge=0;
    double overall_ring_count=0;
    double overall_ring_size=0;
    double overall_hydration=0;

    int cluster_charge;
    double percent_clustered_molecules,strong_ratio;
    FILE *fp_in = NULL, *fp_out=NULL, *fp_top=NULL;
    real time_to_start=0.0;
    int print_every_x_confs=1;
    struct t_fileio* fio = NULL;
    static int verbose_level=0;
    static int check_strict_flag=0;
    static int greater_than_flag=0;
    static int strong_connections_flag=0;
    static int threshold=-1;
    static int threshold_flag=0;
    static int probability_flag=0;
    static int XTC_in_flag=0;
    static int ring_flag=0;
    static int periodicBoundary_flag=1;
    const char *ext = ".xtc";
    int xlen = strlen(ext);
    int slen;
    /*------------------------- START: read the arguments-------------------------*/
    int c;
    while(( c = getopt(argc, argv, "f:o:v:t:cgs:hpm:nra")) != -1 )
    {
        switch(c)
        {
            case 'h':
                printf("Usage: ./test [OPTIONS].. [ARGUMENTS]..\n");
                printf("  -h\t\t: Prints the help information and exits\n");
                printf("  -f <file>\t: Specify the input file\n");
                printf("  -t <file>\t: Specify the .top file in case of XTC input\n");
                printf("  -o <file>\t: Specify the output file\n");
                printf("  -v <int>\t: Specify the verbose level among {1,2,3,4}.\n");
                printf("  -c \t\t: Enables checking if stricter connectedness conditions affects number of connections in each configuration\n");
                printf("  -s <int>\t: Specify size of clusters to be outputted in PDB format\n");
                printf("  -g \t\t: Used with -s to include clusters having size greater than or equal to the argument for -s\n");
                printf("  -p \t\t: Prints the percentage of molecule belonging to a cluster rather than the number of molecules\n");
                printf("  -m <int>\t: Prints statisitics every argument number of configurations(default=1). Needs verbose flag to be disabled\n");
                printf("  -n \t\t: Uses only strong bonds to perform clustering\n");
                printf("  -r \t\t: Performs ring analysis and outputs rings instead\n");
                printf("  -a \t\t: Specifies that the input file is that of a supercell and therefore Periodic Boundary Condition will not applied\n");
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
            case 'n':
                strong_connections_flag=1;
                break;
            case 'r':
                ring_flag=1;
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
            case 'a':
                periodicBoundary_flag=0;
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
    FILE *fp_ring=fopen("ring_statistics.dat","w");
    FILE *fp_cms=fopen("cluster_max_size.dat","w");
    FILE *fp_Kstats=fopen("Kstatistics.dat","w");
    //FILE *fp_testout=fopen("check.pdb","w");
    FILE *fp_csv=fopen("Results.csv","w");
    /*------------------------- END: read the arguments-------------------------*/
    int start_mol_no=-1;
    int no_of_molecules=0;
    int no_of_SOL=0;
    int i,j;

    static ION molecules[MAX_M];
    static COUNTERION Kmolecules[MAX_M];
    static SOL SOLmolecules[MAX_M];
    static int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES];
    static int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES];
    static int SOL_ION_adjacency_matrix[MAX_SOL][MAX_MOLECULES];
    static int cluster_COUNTERION_matrix[MAX_MOLECULES][MAX_MOLECULES];
    static int cluster_SOL_matrix[MAX_MOLECULES][MAX_SOL];

    //Ring analysis
    static int D[MAX_MOLECULES][MAX_MOLECULES];
    static pathArray *P[MAX_MOLECULES][MAX_MOLECULES];
    static pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES];
    vector<ringCandidate> CSet; 
    pathArray CSSSR;


    coordinates boxlength;//, coordinate;
    int conf_number,conf,number_of_K_molecules_in_cluster;
    ION* mol_start=NULL;
    COUNTERION* Kmol_start=NULL;
    SOL* SOLmol_start=NULL;
    /*-------------------------START: read the file --------------------------*/
    if(XTC_in_flag==1)
    {
        if(fp_top==NULL)
        {
            printf("Top file not specified\n");
            exit(0);
        }
        XTC_reader(fio,fp_top,molecules,Kmolecules,SOLmolecules,boxlength,&no_of_molecules,&start_mol_no,&conf_number,time_to_start,&no_of_SOL);

    }
    else
        PDB_reader(fp_in,molecules,Kmolecules,SOLmolecules, boxlength,&no_of_molecules,&start_mol_no,&conf_number, &no_of_SOL);

    int conf_true=conf_number;
    if(verbose_level>=1)
    {
        printf("Number of configurations: %d\n\n",conf_number);
        printf("Number of molecules: %d\n\n",no_of_molecules);
    }
    if (verbose_level>=3)
        printf("BoxLength = %lf %lf %lf \n\n",boxlength[0],boxlength[1],boxlength[2]);

    //printf("no_of_SOL:%d\n",no_of_SOL);
    /*-------------------------END: read the file --------------------------*/
    //Constructing the graph by the means of an adjacency list
    //printf("hi1\n");
    //vector<t_cluster> clusters(no_of_molecules,t_cluster());
    t_cluster clusters[no_of_molecules];
    for(i=0;i<no_of_molecules;i++)
    {
        clusters[i].ION_list=NULL;
        clusters[i].COUTERION_list=NULL;
        clusters[i].SOL_list=NULL;
        clusters[i].ION_list_size=0;
        clusters[i].COUNTERION_list_size=0;
        clusters[i].SOL_list_size=0;
    }
    // for(i=0;i<no_of_molecules;i++)
    //     {
    //         clusters[i].ION_list.reserve(no_of_molecules);
    //         clusters[i].COUTERION_list.reserve(no_of_molecules);
    //         clusters[i].SOL_list.reserve(no_of_molecules);
    //         clusters[i].ringElements.clear();
    //     }
    //printf("hi1\n");
    stack adjacency_list[no_of_molecules];
    int coordination_no[no_of_molecules];
    double average_cluster_coordination_no[no_of_molecules];
    double coordination_no_for_cluster_size[no_of_molecules+1];
    int connectedness[MAX_CONNECTIONS+1];
    int visited[no_of_molecules];
    int number_of_clusters;
    //stack cluster[no_of_molecules];
    int cluster_max_size=0;
    int cluster_size[no_of_molecules+1];
    double cluster_coordination_statistic[no_of_molecules+1];

    vector<ringElements> CSSSR_Elements;
    int number_of_rings;
    int ring_size[no_of_molecules+1];
    int ring_max_size=0;



    //Printing the start of the PDB file
    if(fp_out!=NULL)
    {
        fprintf(fp_out, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
                boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);
        // fprintf(fp_testout, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
        //         boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);
        // printf("CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
        //         boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);     
    }

    //conf_number=1;

    fprintf(fp_csv,"Conf,CNo.,ION,CION,SOL\n");
    
    for(conf=0;conf<conf_true;conf++)
    {
        if (verbose_level>=1)
            printf("Configuration number : %d\n\n",conf+1);
        mol_start=&molecules[conf*no_of_molecules];
        Kmol_start=&Kmolecules[conf*no_of_molecules];
        SOLmol_start=&SOLmolecules[conf*no_of_SOL];
        
        
        //Check if strictness matters
        if(check_strict_flag)
        {
            strict_vs_relaxed(mol_start,boxlength,no_of_molecules, periodicBoundary_flag);
            printf("\n");
        }
        

        if (verbose_level>=4)
        {
            //Prints all the molecules
            
            for(i=0;i<no_of_molecules;i++)
                print_ION(&mol_start[i]);
            printf("\n");
        }


        /*------------------------START: connectedness calculations-----------------*/
        
        /*This part is parallelized */
        //Useful when you dont want to claculate the strong ratio
        //adjacency_complete(mol_start,boxlength,no_of_molecules,adjacency_list,(verbose_level>=3),strong_connections_flag);
        adjacency_matrix_populator(mol_start,boxlength,no_of_molecules,adjacency_matrix,periodicBoundary_flag);
        counterion_adjacency_matrix_populator(mol_start,Kmol_start,boxlength,no_of_molecules,Kadjacency_matrix,periodicBoundary_flag);
        SOL_adjacency_matrix_populator(mol_start,SOLmol_start, boxlength, no_of_molecules, no_of_SOL, SOL_ION_adjacency_matrix, periodicBoundary_flag);
        adjacency_list_from_matrix(adjacency_matrix,no_of_molecules,adjacency_list,(verbose_level>=3),strong_connections_flag);

    // for(j=0;j<no_of_molecules;j++)
    // {
    //     int count=0;
    //     for(i=0;i<no_of_SOL;i++)
    //         if(SOL_ION_adjacency_matrix[i][j])
    //             count++;
    //     if(count>10000)
    //         printf("MOL: %d, SOL:%d\n", j,count);

    // }
        count_counterion_affinity(fp_Kstats, Kadjacency_matrix, no_of_molecules);

        // int countit=0;
        // for(i=0;i<no_of_SOL;i++)
        //     for(j=0;j<no_of_molecules;j++)
        //         if(SOL_ION_adjacency_matrix[i][j])
        //         {
        //             countit++;
        //             //break;
        //         }
        // printf("SOL in conf: %d\n", countit);

        //Ring Analysis
        // if(ring_flag)
        // {
        //     makePIDmatrix(adjacency_matrix, no_of_molecules, D, P, P_dash,(verbose_level>=0));
        //     ringCandidateSearch(&CSet, no_of_molecules, D,P,P_dash);
        //     findSSSR(&CSSSR,&CSet);
        //     // for(auto x: CSet)
        //     //     printf("CNum:%5d | P:%5ld | P_Dash:%5ld\n",x.CNum,x.P->size(),x.P_dash->size());
        //     ringElements* temp;
        //     for(auto x: CSSSR)
        //     {
        //         temp=ringToElements(&x);
        //         CSSSR_Elements.push_back(*temp);
        //         //printRingElements(temp);
        //         delete temp;
        //     }
        //     if(verbose_level>=1)
        //         printRingElementsArray(&CSSSR_Elements);
        //     if(verbose_level>=2)
        //         printPathArray(&CSSSR);

        //     number_of_rings=CSSSR_Elements.size();
        //     for(i=0;i<=no_of_molecules;i++)
        //         ring_size[i]=0;
        //     for(i=0;i<number_of_rings;i++)
        //         ring_size[CSSSR_Elements[i].size()]++;
        //     //smallest ring size is 3
        //     for(i=3;i<=no_of_molecules;i++)
        //         fprintf(fp_ring,"%d ",ring_size[i]);
        //     fprintf(fp_ring,"\n");
        //     overall_ring_count+=number_of_rings;
        //     if(number_of_rings!=0)
        //         ring_max_size=CSSSR_Elements[CSSSR_Elements.size()-1].size();
        //     else
        //         ring_max_size=0;
        //     overall_ring_size+=ring_max_size;

        //     //CSet.clear();
        //     //CSSSR.clear();
        // }


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
            connectedness[i]=0;    

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


        //TODO: Clusters
        //vector<t_cluster>().swap(clusters);
        //printf("%d\n", clusters.size());
        //clusters.clear();
        //clusters.shrink_to_fit();
        //vector<t_cluster> temp;
        //temp.swap(clusters);
        //printf("maybe\n");

        //printf("should be good here\n");

        //#pragma omp parallel for
        for(i=0;i<number_of_clusters;i++)
        {
            //printf("conf:%d\n",conf);
            //printf("ION: %p\n",clusters[i].SOL_list);
            // if(clusters[i].ION_list!=NULL)
            // {
            //     delete[] clusters[i].ION_list;
            //     clusters[i].ION_list=NULL;
            //     //printf("ION dealloc\n");
            // }
            clusters[i].ION_list_size=0;
            // if(clusters[i].COUTERION_list!=NULL)
            // {
            //     delete[] clusters[i].COUTERION_list;
            //     clusters[i].COUTERION_list=NULL;
            //     //printf("CION dealloc\n");
            // }
            clusters[i].COUNTERION_list_size=0;
            // if(clusters[i].SOL_list!=NULL)
            // {
            //     delete[] clusters[i].SOL_list;
            //     clusters[i].SOL_list=NULL;
            //     printf("SOl dealloc\n");
            // }
            clusters[i].SOL_list_size=0;
            clusters[i].ringElements.clear();
        }
        
        
        int cluster_ION_frequency[number_of_clusters]={0};
        for(i=0;i<no_of_molecules;i++)
        {
            cluster_ION_frequency[visited[i]]++;
        }
        // for(i=0;i<number_of_clusters;i++)
        // {
        //     printf("%d ",cluster_ION_frequency[i]);
        // }
        for(i=0;i<number_of_clusters;i++)
        {
            if(clusters[i].ION_list_size<cluster_ION_frequency[i])
            {
                free(clusters[i].ION_list);
                clusters[i].ION_list= (int*) malloc(sizeof(int)*cluster_ION_frequency[i]);
            }
            //clusters[i].ION_list= (int*)malloc(sizeof(int)*cluster_ION_frequency[i]);
        }
        //printf("should be bad here\n");
            
        //Adding molecule number to ION_list
        for(i=0;i<no_of_molecules;i++)
        {
            clusters[visited[i]].ION_list[clusters[visited[i]].ION_list_size++]=i;

        }
        //printf("should be bad here\n");

        // for(i=0;i<number_of_clusters;i++)
        // {
        //     for (j=0;j<clusters[i].ION_list_size;j++)
        //         printf("%d ",clusters[i].ION_list[j]);
        //     printf("\n");
        // }
        //printf("maybe\n");
        //printf("hi1\n");

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
            

        // //Initializing the stacks
        // for(i=0;i<number_of_clusters;i++)
        // {
        //     cluster[i].length=0;
        //     cluster[i].top=NULL;
        // }

        // //Pushing molecule number on the stack corresponding to the cluster number
        // for(i=0;i<no_of_molecules;i++)
        //     add_node_given_value(&cluster[visited[i]],i);

        

        //Recording cluster_size and cluster_max_size
        cluster_max_size=0;
        // for(i=0;i<number_of_clusters;i++)
        // {
        //     cluster_size[cluster[i].length]++;
        //     if(cluster[i].length>cluster_max_size)
        //         cluster_max_size=cluster[i].length;
        // }
        //TODO: CLusters
        for(i=0;i<number_of_clusters;i++)
        {
            cluster_size[clusters[i].ION_list_size]++;
            if(clusters[i].ION_list_size>cluster_max_size)
                cluster_max_size=clusters[i].ION_list_size;
        }
        //printf("hi1\n");

        /*
        //TODO: Obsoltet code

        //Calculation average coordination number of cluster
       // calculate_cluster_coordination_number(cluster,coordination_no,average_cluster_coordination_no,number_of_clusters);

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


        //TODO: Obsoltet code END
        */

        //printf("hi1\n");

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
        //printf("hi1\n");

        add_COUNTERION_to_cluster(Kadjacency_matrix, cluster_COUNTERION_matrix,no_of_molecules, no_of_molecules,clusters,number_of_clusters);
        //printf("Done with cion\n");


        add_SOL_to_cluster(SOL_ION_adjacency_matrix, cluster_SOL_matrix,no_of_molecules, no_of_SOL,clusters,number_of_clusters);
        // int count=0;
        // for(i=0;i<no_of_molecules;i++)
        //     for(j=0;j<no_of_molecules;j++)
        //         if(Kadjacency_matrix[i][j])
        //         {
        //             count+=1;
        //             break;
        //         }
        // printf("Acutal: %d\n",count);


        if(ring_flag)
        {
            //printf("Conf for ring:%d\n",conf);
            number_of_rings=0;
            ring_max_size=0;
            for(i=0;i<=no_of_molecules;i++)
                ring_size[i]=0;
            for(i=0;i<number_of_clusters;i++)
            {
                if(clusters[i].ION_list_size>=3)
                {
                    ringDriver(adjacency_matrix,clusters[i].ION_list,clusters[i].ION_list_size,D,P,P_dash,strong_connections_flag, verbose_level,&CSet,&CSSSR,&(clusters[i].ringElements));
                    //printf("Number of rings:%d\n",clusters[i].ringElements.size());
                    number_of_rings+=clusters[i].ringElements.size();
                    for(auto x: clusters[i].ringElements)
                    {
                        ring_size[x.size()]++;
                        if(x.size()>ring_max_size)
                            ring_max_size=x.size();

                        //     //smallest ring size is 3
                        //     for(i=3;i<=no_of_molecules;i++)
                        //         fprintf(fp_ring,"%d ",ring_size[i]);
                        //     fprintf(fp_ring,"\n");
                        //     overall_ring_count+=number_of_rings;
                        //     if(number_of_rings!=0)
                        //         ring_max_size=CSSSR_Elements[CSSSR_Elements.size()-1].size();
                        //     else
                        //         ring_max_size=0;
                        //     overall_ring_size+=ring_max_size;
                    }
                }
            }
            for(i=3;i<=no_of_molecules;i++)
                fprintf(fp_ring,"%d ",ring_size[i]);
            fprintf(fp_ring,"\n");

            overall_ring_size+=ring_max_size;
            overall_ring_count+=number_of_rings;
        }
        /*--------------------------END: clustering calculations------------------*/

        /*-----------------------START: Output PDB------------------*/
        for(i=0;i<number_of_clusters;i++)
            fprintf(fp_csv,"%d, %d, %d, %d, %d\n",conf,i,clusters[i].ION_list_size,clusters[i].COUNTERION_list_size,clusters[i].SOL_list_size);
        
        if(threshold_flag==0)
            threshold=cluster_max_size;

        // for(i=0;i<number_of_clusters;i++)
        // {
        //     if(((greater_than_flag)&&(clusters[i].ION_list_size>=threshold))||((!greater_than_flag)&&(clusters[i].ION_list_size==threshold)))
        //         printf("Cluster %d: ION: %d CION: %d SOL: %d\n", i, clusters[i].ION_list_size, clusters[i].COUNTERION_list_size, clusters[i].SOL_list_size);

        // }
        if(fp_out!=NULL)
        {
            //number_of_K_molecules_in_cluster=0;
            //number_of_K_molecules_in_cluster=fprintf_K_ions_in_cluster(fp_out,Kadjacency_matrix,Kmol_start,cluster,no_of_molecules,number_of_clusters,threshold,greater_than_flag);
            //printf("No of K mols in cluster =%d\n",number_of_K_molecules_in_cluster);
            //fprintf_conf_PDB(fp_out,mol_start,cluster,number_of_clusters,threshold,greater_than_flag,number_of_K_molecules_in_cluster);
            fprintf_all(fp_out,mol_start,Kmol_start, SOLmol_start, clusters, number_of_clusters, threshold, greater_than_flag);

        }



        int number_of_K_molecules_in_cluster=0;
        int hydration=0;
        for(i=0;i<number_of_clusters;i++)
            if((greater_than_flag&&clusters[i].ION_list_size>=threshold)||(!greater_than_flag&&clusters[i].ION_list_size==threshold))
            {
                number_of_K_molecules_in_cluster+=clusters[i].COUNTERION_list_size;
                hydration+=clusters[i].SOL_list_size;
            }
                
        //printf("%d %d %d\n",number_of_K_molecules_in_cluster,CION_sum, cluster_size[threshold]);
        //assert(CION_sum==number_of_K_molecules_in_cluster);


        /*-----------------------END: Output PDB------------------*/

        /*-----------------------START: Print Stats------------------*/
        // if(threshold_flag==0)
        //     threshold=cluster_max_size;
        // if(fp_out!=NULL)
        // {
        //     fprintf_conf_PDB(fp_out,mol_start,cluster,number_of_clusters,threshold,greater_than_flag);
        // }

        /*-----------------------END: Print Stats------------------*/


        overall_cluster_size+=cluster_max_size;
        percent_clustered_molecules= ((double)(no_of_molecules-cluster_size[1]))/no_of_molecules*100;
        overall_percentage_clustered+=percent_clustered_molecules;
        strong_ratio=strong_connection_ratio(adjacency_matrix,no_of_molecules)*100;
        overall_percentage_strong+=strong_ratio;
        overall_hydration+=hydration;
        
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
            if(ring_flag)
                printf("Conf: %5d | MaxClusterSize: %5d | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %5d | Hydration : %5d | MaxRingSize : %5d | RingCount : %5d\n",conf,cluster_max_size,percent_clustered_molecules,strong_ratio,cluster_charge,hydration,ring_max_size,number_of_rings);
            else
                printf("Conf: %5d | MaxClusterSize: %5d | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %5d | Hydration : %5d\n",conf,cluster_max_size,percent_clustered_molecules,strong_ratio,cluster_charge, hydration);
        }

        /*-----------------------START: Cleanup------------------*/
        // if(ring_flag)
        // {
        //     printf("Conf for ring:%d\n",conf);
        //     for(i=0;i<no_of_molecules;i++)
        //         for(j=0;j<no_of_molecules;j++)
        //         {
        //             delete P[i][j];
        //             delete P_dash[i][j];
        //         }
        //     printf("Conf for ring:%d\n",conf);
        // }



        // for(i=0;i<number_of_clusters;i++)
        // {
        //     empty_stack(&cluster[i]);
        // }
        
        CSSSR_Elements.clear();
        CSet.clear();
        CSSSR.clear();
        // for(i=0;i<number_of_clusters;i++)
        // {
        //     //TODO: Can further optimise by just reallocating when size is less

        //     //printf("conf:%d\n",conf);
        //     //printf("ION: %p\n",clusters[i].SOL_list);


        //     if(clusters[i].ION_list!=NULL)
        //     {
        //         free(clusters[i].ION_list);
        //         clusters[i].ION_list=NULL;
        //         //printf("ION dealloc\n");
        //     }
        //     clusters[i].ION_list_size=0;
        //     if(clusters[i].COUTERION_list!=NULL)
        //     {
        //         free(clusters[i].COUTERION_list);
        //         clusters[i].COUTERION_list=NULL;
        //         //printf("CION dealloc\n");
        //     }
        //     clusters[i].COUNTERION_list_size=0;
        //     if(clusters[i].SOL_list!=NULL)
        //     {
        //         free(clusters[i].SOL_list);
        //         clusters[i].SOL_list=NULL;
        //         //printf("SOl dealloc\n");
        //     }
        //     clusters[i].SOL_list_size=0;
        //     clusters[i].ringElements.clear();
        // }



        // printf("Deleting teh lists\n");
        // for(i=0;i<number_of_clusters;i++)
        // {
        //     delete clusters[i].ION_list;
        //     printf("Deleted ion teh lists\n");
        //     delete clusters[i].COUTERION_list;
        //     printf("Deleted cion lists\n");
        //     delete clusters[i].SOL_list;
        //     printf("Deleted sol lists\n");
        // }
        //break;

        /*-----------------------END: Cleanup------------------*/

    //printf("hi %d %d \n",conf_number, conf_true);
    }
    printf("Exited loop\n");
    overall_cluster_size/=conf_true;
    overall_percentage_clustered/=conf_true;
    overall_percentage_strong/=conf_true;
    overall_cluster_charge/=conf_true;
    overall_ring_count/=conf_true;
    overall_ring_size/=conf_true;
    overall_hydration/=conf_true;

    if(!ring_flag)
        printf("\nOverall     | MaxClusterSize: %5.2lf | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %5.2lf | Hydration : %5.2lf \n",overall_cluster_size,overall_percentage_clustered,overall_percentage_strong,overall_cluster_charge, overall_hydration);
    else
        printf("\nOverall     | MaxClusterSize: %5.2lf | %%age Clustered: %5.2lf | %%age Strong : %5.2lf | Cluster Charge : %5.2lf | Hydration : %5.2lf | MaxRingSize : %5.2lf | RingCount : %5.2lf\n",overall_cluster_size,overall_percentage_clustered,overall_percentage_strong,overall_cluster_charge,overall_hydration,overall_ring_size,overall_ring_count);
    if(fp_out!=NULL)
        fclose(fp_out);

    fclose(fp_stats);
    fclose(fp_cms);
    fclose(fp_ring);
    fclose(fp_Kstats);
    //fclose(fp_testout);
    fclose(fp_csv);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
 
    printf("The elapsed time is %f seconds\n", elapsed);
    return 0;
}