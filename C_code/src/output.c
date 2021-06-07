#include "output.h"

void fprintf_cluster_PDB(FILE* fp_out,HPO molecules[],stack* cluster,int *mol_no,int* atom_no)
{
    node *temp=cluster->top;
    while(temp!=NULL)
    {
        (*mol_no)++;
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  PL  HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].P[0],molecules[temp->data].P[1],molecules[temp->data].P[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  OHL HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].OHL_1[0],molecules[temp->data].OHL_1[1],molecules[temp->data].OHL_1[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  HOL HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].HOL_1[0],molecules[temp->data].HOL_1[1],molecules[temp->data].HOL_1[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  OHL HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].OHL_2[0],molecules[temp->data].OHL_2[1],molecules[temp->data].OHL_2[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  HOL HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].HOL_2[0],molecules[temp->data].HOL_2[1],molecules[temp->data].HOL_2[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  O2L HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].O2L_1[0],molecules[temp->data].O2L_1[1],molecules[temp->data].O2L_1[2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  O2L HPO X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), (*mol_no),molecules[temp->data].O2L_2[0],molecules[temp->data].O2L_2[1],molecules[temp->data].O2L_2[2]);
        temp=temp->next;
    }

}

void fprintf_conf_PDB(FILE* fp_out,HPO molecules[],stack clusters[],int number_of_clusters, int threshold, int greater_than_flag)
{
        int mol_no=0;
        int atom_no=0;
        int i;
        if(greater_than_flag==1)
        {
                for(i=0; i<number_of_clusters; i++)
                        if(clusters[i].length>=threshold)
                                fprintf_cluster_PDB(fp_out,molecules,&clusters[i],&mol_no,&atom_no);
        }
        else
                for(i=0; i<number_of_clusters; i++)
                        if(clusters[i].length==threshold)
                                fprintf_cluster_PDB(fp_out,molecules,&clusters[i],&mol_no,&atom_no);
        fprintf(fp_out, "END\n");
}