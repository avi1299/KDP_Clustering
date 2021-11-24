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
