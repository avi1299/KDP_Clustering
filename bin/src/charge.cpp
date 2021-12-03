#include "charge.h"
using namespace std;

//int parallelism_enabled=1;

int fprintf_K_ions_in_cluster(FILE* fp_out, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , COUNTERION Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag){
    int i,j;
    int K_mols_of_interest[no_of_molecules]={0};

    node * temp=NULL;
    for(i=0;i<number_of_clusters;i++)
    {
        if(((greater_than_flag)&&(clusters[i].length>=threshold))||((!greater_than_flag)&&(clusters[i].length==threshold)))
        {
            temp=clusters[i].top;
            while(temp!=NULL)
            {
                for(j=0;j<no_of_molecules;j++)
                {
                    K_mols_of_interest[j]=K_mols_of_interest[j]||Kadjacency_matrix[j][temp->data];
                }
                temp=temp->next;
            }
        }
    }

    j=1;
    //Printing the COUNTERION molecules
    int count;
    for(i=0;i<no_of_molecules;i++)
    {
        if(K_mols_of_interest[i])
        {
            count=0;
            for(int k=0;k<no_of_molecules;k++)
                count+=Kadjacency_matrix[i][k];
            fprintf(fp_out, "ATOM  %5d  K   K%d  X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (j),(count), (j),Kmolecules[i].posn[0],Kmolecules[i].posn[1],Kmolecules[i].posn[2]);
            j++;
        }
    }
    j--;
    return j;
            
}

void set_mols_of_interest(int mols_of_interest[], stack *cluster){
    node* temp=cluster->top;
    while(temp!=NULL)
    {
        mols_of_interest[temp->data]=1;
        temp=temp->next;
    }
}

void count_counterion_affinity(FILE* fp_Kstat, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int no_of_molecules){
    vector<int> count(no_of_molecules+1,0);
    #pragma omp parallel for
    for(int i=0;i<no_of_molecules;i++)//Looping over each K
    {
        for(int j=0;j<no_of_molecules;j++)
        {
            count[i]+=Kadjacency_matrix[i][j];
        }
    }
    vector<int> stats(no_of_molecules+1,0);
    for(int i=0;i<no_of_molecules+1;i++)
        stats[count[i]]++;
    for(int i=0;i<no_of_molecules+1;i++)
        fprintf(fp_Kstat,"%d ",stats[i]);
    fprintf(fp_Kstat,"\n");
}

int fprintf_SOL(FILE* fp_out, int SOLadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , SOL SOLmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int no_of_SOL, int current_atom, int current_mol,int threshold, int greater_than_flag){
    int i,j;
    int SOL_mols_of_interest[no_of_SOL]={0};

    node * temp=NULL;
    for(i=0;i<number_of_clusters;i++)
    {
        if(((greater_than_flag)&&(clusters[i].length>=threshold))||((!greater_than_flag)&&(clusters[i].length==threshold)))
        {
            temp=clusters[i].top;
            while(temp!=NULL)
            {
                for(j=0;j<no_of_SOL;j++)
                {
                    SOL_mols_of_interest[j]=SOL_mols_of_interest[j]||SOLadjacency_matrix[j][temp->data];
                }
                temp=temp->next;
            }
        }
    }

    //printf("hehe\n");
    j=current_mol+1;
    //Printing the COUNTERION molecules
    int count;
    int s=current_atom+1;
    //printf("%d\n",no_of_SOL);
    for(i=0;i<no_of_SOL;i++)
    {
        //printf("hi");
        if(SOL_mols_of_interest[i])
        {
            count=0;
            for(int k=0;k<no_of_molecules;k++)
                count+=SOLadjacency_matrix[i][k];
            fprintf(fp_out, "ATOM  %5d  OW  SOL X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (s++)/*,(count)*/, (j),SOLmolecules[i].posn[OW][0],SOLmolecules[i].posn[OW][1],SOLmolecules[i].posn[OW][2]);
            fprintf(fp_out, "ATOM  %5d  HW1 SOL X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (s++)/*,(count)*/, (j),SOLmolecules[i].posn[HW1][0],SOLmolecules[i].posn[HW1][1],SOLmolecules[i].posn[HW1][2]);
            fprintf(fp_out, "ATOM  %5d  HW2 SOL X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (s++)/*,(count)*/, (j),SOLmolecules[i].posn[HW2][0],SOLmolecules[i].posn[HW2][1],SOLmolecules[i].posn[HW2][2]);                
            j++;
            //printf("%d",j);
        }
    }
    j--;
    fprintf(fp_out, "END\n");
    return j;
}

void add_COUNTERION_to_cluster(int CIONadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int cluster_COUNTERION_matrix[MAX_MOLECULES][MAX_MOLECULES], int no_of_ION, int no_of_CION, t_cluster* clusters, int number_of_clusters)
{
    
    int clusterList_size=number_of_clusters;
    int i,j,k;

    //clearing the matrix
    #pragma omp parallel for private(j)
    for(i=0;i<clusterList_size;i++)
    {
        //(*clusters)[i].COUTERION_list.clear();
        for(j=0;j<no_of_CION;j++)
            cluster_COUNTERION_matrix[i][j]=0;
    }
    //memset(cluster_COUNTERION_matrix, 0, sizeof(MAX_MOLECULES*MAX_MOLECULES)); // for automatically-allocated arrays


    //Populating Matrix
    #pragma omp parallel for private(j,k)
    for(i=0;i<clusterList_size;i++)
        for(j=0;j<no_of_CION;j++)
            for(k=0;k<clusters[i].ION_list_size;k++)
                cluster_COUNTERION_matrix[i][j]+=CIONadjacency_matrix[j][clusters[i].ION_list[k]];
    
    int maxnum[no_of_CION];
    int COUNTERION_belongs_to[no_of_CION];
    
    for(i=0;i<no_of_CION;i++)
    {
        maxnum[i]=0;
        COUNTERION_belongs_to[i]=-1;
    }

    #pragma omp parallel for private(i)
    for(j=0;j<no_of_CION;j++)
        for(i=0;i<clusterList_size;i++)
            if(cluster_COUNTERION_matrix[i][j]>maxnum[j])
            {
                maxnum[j]=cluster_COUNTERION_matrix[i][j];
                COUNTERION_belongs_to[j]=i;
            }

    int cluster_COUNTERION_frequency[number_of_clusters]={0};
    for(i=0;i<no_of_CION;i++)
    {
        if(COUNTERION_belongs_to[i]>=0)
            cluster_COUNTERION_frequency[COUNTERION_belongs_to[i]]++;
    }
    
    
    // for(i=0;i<no_of_CION;i++)
    //     printf("%d ",COUNTERION_belongs_to[i]);
    for(i=0;i<number_of_clusters;i++)
    {
        //if(cluster_COUNTERION_frequency[i]>0)
        if(clusters[i].COUNTERION_list_size<cluster_COUNTERION_frequency[i])
            {
                free(clusters[i].COUTERION_list);
                clusters[i].COUTERION_list= (int*) malloc(sizeof(int)*cluster_COUNTERION_frequency[i]);
            }
        //lusters[i].COUTERION_list= (int*) malloc(sizeof(int)*cluster_COUNTERION_frequency[i]);
        clusters[i].COUNTERION_list_size=0;
        // else 
        //     clusters[i].COUTERION_list=new int[1];
        //printf("CION Pointer: %p + %d\n",clusters[i].SOL_list,cluster_COUNTERION_frequency[i]);
    }
    



    for(i=0;i<no_of_CION;i++)
    {
        //printf("%d %d\n",i, COUNTERION_belongs_to[i]);
        if(COUNTERION_belongs_to[i]>=0)
        {
            //printf("%d\n",i);
            //int *index=&(clusters[COUNTERION_belongs_to[i]].COUNTERION_list_size);
            //printf("%d\n",*index);
            clusters[COUNTERION_belongs_to[i]].COUTERION_list[clusters[COUNTERION_belongs_to[i]].COUNTERION_list_size]=i;
            //printf("done\n");
            clusters[COUNTERION_belongs_to[i]].COUNTERION_list_size++;

        }
    }
    //printf("why\n");

    // for(i=0;i<number_of_clusters;i++)
    // {
    //     for (j=0;j<clusters[i].COUNTERION_list_size;j++)
    //         printf("%d ",clusters[i].COUTERION_list[j]);
    //     printf("\n");
    // }
    //printf("here\n");
}

void add_SOL_to_cluster(int SOL_ION_adjacency_matrix[MAX_MOLECULES][MAX_MOLECULES], int cluster_SOL_matrix[MAX_MOLECULES][MAX_MOLECULES], int no_of_ION, int no_of_SOL, t_cluster* clusters, int number_of_clusters)
{
    int clusterList_size=number_of_clusters;
    int i,j,k;

    int count=0;

    #pragma omp parallel for private(j)
    for(i=0;i<clusterList_size;i++)
    {
        //(*clusters)[i].SOL_list.clear();
        for(j=0;j<no_of_SOL;j++)
            cluster_SOL_matrix[i][j]=0;
    }


    #pragma omp parallel for private(j,k)
    for(i=0;i<clusterList_size;i++)
        for(j=0;j<no_of_SOL;j++)
            for(k=0;k<clusters[i].ION_list_size;k++)
                cluster_SOL_matrix[i][j]+=SOL_ION_adjacency_matrix[j][clusters[i].ION_list[k]];
    
    for(i=0;i<clusterList_size;i++)
        for(j=0;j<no_of_SOL;j++)
            count+=cluster_SOL_matrix[i][j];
    //printf("count: %d\n",count);
    int maxnum[no_of_SOL]={0};
    int SOL_belongs_to[no_of_SOL];

    for(i=0;i<no_of_SOL;i++)
    {
        //maxnum[i]=0;
        SOL_belongs_to[i]=-1;
    }

    #pragma omp parallel for private(i)
    for(j=0;j<no_of_SOL;j++)
        for(i=0;i<clusterList_size;i++)
            if(cluster_SOL_matrix[i][j]>maxnum[j])
            {
                maxnum[j]=cluster_SOL_matrix[i][j];
                SOL_belongs_to[j]=i;
            }

    int cluster_SOL_frequency[number_of_clusters]={0};
    for(i=0;i<no_of_SOL;i++)
    {
        if(SOL_belongs_to[i]>=0)
        cluster_SOL_frequency[SOL_belongs_to[i]]++;
    }
    //for(i=0;i<number_of_clusters;i++)
        //printf("space:%d\n", cluster_SOL_frequency[i]);

    for(i=0;i<number_of_clusters;i++)
    {
        
        if(cluster_SOL_frequency[i]>0)
            if(clusters[i].SOL_list_size<cluster_SOL_frequency[i])
            {
                free(clusters[i].SOL_list);
                clusters[i].SOL_list= (int*) malloc(sizeof(int)*cluster_SOL_frequency[i]);
            }
            //assert(clusters[i].SOL_list!=NULL);
            clusters[i].SOL_list_size=0;
        // else
        //     clusters[i].SOL_list=new int[1];
        //printf("SOL Pointer: %p + %d\n",clusters[i].SOL_list,cluster_SOL_frequency[i]);
    }



    //printf("it should be wrong here\n");
    for(i=0;i<no_of_SOL;i++)
    {
        //printf("%d %d\n",i, COUNTERION_belongs_to[i]);
        if(SOL_belongs_to[i]!=-1)
        {
            //(*clusters)[SOL_belongs_to[i]].SOL_list.push_back(i);
            //int *index=&(clusters[SOL_belongs_to[i]].SOL_list_size);
            //printf("%d\n",*index);
            clusters[SOL_belongs_to[i]].SOL_list[clusters[SOL_belongs_to[i]].SOL_list_size]=i;
            //printf("%d\n",*index);
            clusters[SOL_belongs_to[i]].SOL_list_size++;
        }
    }
    // for(i=0;i<number_of_clusters;i++)
    // {
    //     for (j=0;j<clusters[i].SOL_list_size;j++)
    //         printf("%d ",clusters[i].SOL_list[j]);
    //     printf("\n");
    // }
    //printf("this ran\n");
}