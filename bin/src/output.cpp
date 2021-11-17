#include "output.h"

void fprintf_cluster_PDB(FILE* fp_out,ION molecules[],stack* cluster,int *mol_no,int* atom_no, int index)
{
    node *temp=cluster->top;
    while(temp!=NULL)
    {
        (*mol_no)++;
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  PL  HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[PL][0],molecules[temp->data].posn[PL][1],molecules[temp->data].posn[PL][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  OHL HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[OHL_1][0],molecules[temp->data].posn[OHL_1][1],molecules[temp->data].posn[OHL_1][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  HOL HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[HOL_1][0],molecules[temp->data].posn[HOL_1][1],molecules[temp->data].posn[HOL_1][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  OHL HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[OHL_2][0],molecules[temp->data].posn[OHL_2][1],molecules[temp->data].posn[OHL_2][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  HOL HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[HOL_2][0],molecules[temp->data].posn[HOL_2][1],molecules[temp->data].posn[HOL_2][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  O2L HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[O2L_1][0],molecules[temp->data].posn[O2L_1][1],molecules[temp->data].posn[O2L_1][2]);
        (*atom_no)++;
        fprintf(fp_out, "ATOM  %5d  O2L HP%d X%4d   %8.3lf%8.3lf%8.3lf  0.00  0.00\n",
                (*atom_no), index, (*mol_no),molecules[temp->data].posn[O2L_2][0],molecules[temp->data].posn[O2L_2][1],molecules[temp->data].posn[O2L_2][2]);
        temp=temp->next;
    }

}

void fprintf_conf_PDB(FILE* fp_out,ION molecules[],stack clusters[],int number_of_clusters, int threshold, int greater_than_flag, int HPO_start_no)
{
        int mol_no=HPO_start_no;
        int atom_no=HPO_start_no;
        int i;
        int index=0;
        if(greater_than_flag==1)
        {
                for(i=0; i<number_of_clusters; i++)
                        if(clusters[i].length>=threshold)
                                fprintf_cluster_PDB(fp_out,molecules,&clusters[i],&mol_no,&atom_no,index++);
        }
        else
                for(i=0; i<number_of_clusters; i++)
                        if(clusters[i].length==threshold)
                                fprintf_cluster_PDB(fp_out,molecules,&clusters[i],&mol_no,&atom_no,index++);
        fprintf(fp_out, "END\n");
}