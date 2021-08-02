#include "charge.h"
using namespace std;

//int parallelism_enabled=1;

int fprintf_K_ions_in_cluster(FILE* fp_out, HPO molecules[], K Kmolecules[], stack clusters[], int no_of_molecules, int number_of_clusters, int threshold, int greater_than_flag, coordinates boxlength){
    int i,j;
    int HPO_mols_of_interest[no_of_molecules]={0},K_mols_of_interest[no_of_molecules]={0};
    if(greater_than_flag==1)
    {
            for(i=0; i<number_of_clusters; i++)
                    if(clusters[i].length>=threshold)
                            set_mols_of_interest(HPO_mols_of_interest,&clusters[i]);
    }
    else
            for(i=0; i<number_of_clusters; i++)
                    if(clusters[i].length==threshold)
                            set_mols_of_interest(HPO_mols_of_interest,&clusters[i]);
    
    //#pragma omp parallel for schedule(static, 1) private(j) if (parallelism_enabled)
    for(i=0;i<no_of_molecules;i++)
    {
        if(HPO_mols_of_interest[i])
        {
            for(j=0;j<no_of_molecules;j++)
            {
                if(connected_K_HPO(&Kmolecules[j],&molecules[i],boxlength))
                    K_mols_of_interest[j]=1;
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

