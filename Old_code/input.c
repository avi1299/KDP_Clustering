/*I assume each of the molecules in the system is in a cluster, then I merge the cluster by its connection*/
/*By lots of iteration, we can find all the seperate cluster*/
/*This code read a coordinate file and outpu the larest cluster */
#include "input.h"
int main(int argc,char *argv[])
{

        FILE *fp_in, *fp_out;
        char filename[NAME], xname[2][NAME], line[LLEN];
        int i, j, m, n, water_nmols=0, mol_tag, mol_type, classify, go_or_not, connect_or_not, standard=0;
        double box[3][2]={0.0}, temp;
        static HPO water[MAX_M]; //only for 50000 water molecules
        static cluster_mt cluster[MAX_M];//for all the clusters
        coordinates boxlength, coordinate;
        /*------------------------------ read the arguments-----------------------------------------*/
        int c;
        while(( c = getopt(argc, argv, "f:o:")) != -1 )
        {
                switch(c)
                {
                case 'f':
                        if(( fp_in  =fopen(optarg,"r"))==NULL)
                        {printf("cannot open Infile \n"); exit(0); }
                        break;
                case 'o':
                        if(( fp_out  =fopen(optarg,"w"))==NULL)
                        {printf("cannot open Outfile \n"); exit(0); }
                        break;
                }
        }
/*-------------------------------read the inout file----------------------------------------*/
        if(fgets(line, LLEN, fp_in) == NULL)              //the fist line
                printf("Error: infile is empty!\n");
        fgets(line, LLEN, fp_in);                //the second line
        fgets(line, LLEN, fp_in);                // the third line NUMBER of Atoms
        fgets(line, LLEN, fp_in);                //number of Atoms
        sscanf(line, "%d", &water_nmols);
        fgets(line, LLEN, fp_in);                //BOX BONUS
        for(i=0; i<3; i++)
        {
                fgets(line, LLEN, fp_in);       for(;; )
        {
                if((fgets(line, LLEN, fp_in) == NULL))
                        break;
                if(sscanf(line, "%d%d%lf%lf%lf%d", &mol_tag, &mol_type,//mol_tag, mol_type,
                          &coordinate[0], &coordinate[1], &coordinate[2], &classify) != 6) break;
                if(mol_type == 1)//molecule type
                {
                        water[mol_tag-1].tag           =  mol_tag;//tag alway begin with 0
                        for(i=0; i<3; i++)
                                water[mol_tag-1].oxygen[i]   =  coordinate[i]*boxlength[i];
                        cluster[mol_tag-1].nodeid[0] = mol_tag -1;
                        cluster[mol_tag-1].node_number = 0;// initilizing with 0
                        if(classify==1 || classify==2)
                        {
                                cluster[mol_tag-1].node_number = 1;
                                standard++;
                        }
                }
        }
        printf("The system contains %d ice molecules
                sscanf(line, "%lf%lf", &box[i][0], &box[i][1]);
                boxlength[i] = box[i][1] - box[i][0];
        }
        fgets(line, LLEN, fp_in);               //skip the blank line before the coordinates
        for(;; )
        {
                if((fgets(line, LLEN, fp_in) == NULL))
                        break;
                if(sscanf(line, "%d%d%lf%lf%lf%d", &mol_tag, &mol_type,//mol_tag, mol_type,
                          &coordinate[0], &coordinate[1], &coordinate[2], &classify) != 6) break;
                if(mol_type == 1)//molecule type
                {
                        water[mol_tag-1].tag           =  mol_tag;//tag alway begin with 0
                        for(i=0; i<3; i++)
                                water[mol_tag-1].oxygen[i]   =  coordinate[i]*boxlength[i];
                        cluster[mol_tag-1].nodeid[0] = mol_tag -1;
                        cluster[mol_tag-1].node_number = 0;// initilizing with 0
                        if(classify==1 || classify==2)
                        {
                                cluster[mol_tag-1].node_number = 1;
                                standard++;
                        }
                }
        }
        printf("The system contains %d ice molecules\n", standard);
/*-----------------------CLustering the system----------------------------------------------*/
        for(;; )
        {
                go_or_not = 0;
                for(i=0; i<water_nmols; i++)
                        for(j=0; j<water_nmols; j++)
                        {
                                if(i!=j&&cluster[i].node_number>0&&cluster[j].node_number>0)
                                {
                                        connect_or_not=0;
                                        for(m=0; m<cluster[i].node_number; m++)
                                                for(n=0; n<cluster[j].node_number; n++)
                                                {
                                                        temp = mindist(water[cluster[i].nodeid[m]].oxygen, water[cluster[j].nodeid[n]].oxygen, boxlength);
                                                        if(temp<=CUTOFF)
                                                        {connect_or_not =1;
                                                         break;}
                                                }
                                        if(connect_or_not==1)
                                        {
                                                for(m=cluster[i].node_number; m<cluster[i].node_number+cluster[j].node_number; m++)
                                                        cluster[i].nodeid[m] = cluster[j].nodeid[m-cluster[i].node_number];
                                                cluster[i].node_number += cluster[j].node_number;
                                                cluster[j].node_number = 0;
                                                go_or_not = 1;
                                        }
                                }
                        }
                if(go_or_not==0)
                        break;
        }
/*------------------------------------------order the cluster------------------------------------*/
        int largest=0, temp1=0;
        for(i=0; i<water_nmols; i++)
                if(cluster[i].node_number>temp1)
                {
                        largest= i;
                        temp1 = cluster[i].node_number;
                }
/*------------------------------output the pdb file------------------------------------------*/
//   fprintf(fp_out, "TITLE    A system configuration for the Lammps datafile\n");
//   fprintf(fp_out,  "REMARK   THIS IS A SIMULATION BOX\n");
        fprintf(fp_out, "CRYST1%9.3lf%9.3lf%9.3lf%7.2lf%7.2lf%7.2lf P 1           1\n",
                boxlength[0], boxlength[1], boxlength[2], 90.0, 90.0, 90.0);
//   fprintf(fp_out, "MODEL        1\n");
        for(i=0; i<cluster[largest].node_number; i++)
                fprintf(fp_out, "ATOM  %5d   O  SOL  %4d    %8.3lf%8.3lf%8.3lf  1.00  0.00\n",
                        cluster[largest].nodeid[i], cluster[largest].nodeid[i],
                        water[cluster[largest].nodeid[i]].oxygen[0], water[cluster[largest].nodeid[i]].oxygen[1],
                        water[cluster[largest].nodeid[i]].oxygen[2]);
        fprintf(fp_out, "TER\n");
        fprintf(fp_out, "ENDMDL\n");
        fclose(fp_in);
        fclose(fp_out);
        return 0;
}
