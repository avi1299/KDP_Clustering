#include "charge.h"
using namespace std;

//int parallelism_enabled=1;

int fprintf_K_ions_in_cluster(FILE* fp_out, int Kadjacency_matrix[MAX_MOLECULES][MAX_MOLECULES] , K Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag){
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
    for(i=0;i<no_of_molecules;i++)
    {
        if(K_mols_of_interest[i])
        {
            fprintf(fp_out, "ATOM  %5d  K   K   X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (j), (j),Kmolecules[i].posn[0],Kmolecules[i].posn[1],Kmolecules[i].posn[2]);
            j++;
        }
    }
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

