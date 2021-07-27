#include "input.h"
#include "gromacs/fileio/xtcio.h"

//Assumption all the molecules are listed together and though order of molecules doesn't matter, all the 7 atoms corresponding to a molecule should be together
//i.e. one after the other.
void PDB_reader(FILE* fp_in,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number)
{
    int i,atom_no,mol_no,old_mol_no=-1;;
    char line[LLEN],atom_name[LLEN],mol_name[LLEN],tag[LLEN];
    int HOL_read=0;
    int OHL_read=0;
    int O2L_read=0;
    int k_count=0;
    *no_of_molecules=-1;
    *conf_number=0;
    coordinates coordinate;
    if(fgets(line, LLEN, fp_in) == NULL)
    {
        printf("Error: infile is empty!\n"); 
        exit(1);
    }
    sscanf(line,"%*s %lf %lf %lf %*lf %*lf %*lf %*s %*d %*d",&boxlength[0],&boxlength[1],&boxlength[2]);
    //printf("%s",line);
    //printf("BoxLength = %lf %lf %lf \n",boxlength[0],boxlength[1],boxlength[2]);
    for(;;)
    {
        if(fgets(line, LLEN, fp_in) == NULL)
            break;
        sscanf(line,"%s %d %s %s X %d %lf %lf %lf %*lf %*lf",tag,&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
        if(strcmp(tag,"END")==0)
            (*conf_number)++;
        //Executes below when it encounters HPO molecules. It extracts information for every atom in every HPO molecule 
        else if(strcmp(mol_name,"K")==0)
        {
            for(i=0;i<3;i++)
                    Kmolecules[k_count].posn[i]=coordinate[i];
            k_count++;
        }
        else if(strcmp(mol_name,"HPO")==0)
        {
            //Noting where the first HPO molecule was encountered
            if(*start_mol_no==-1)
                *start_mol_no=mol_no;
            if(old_mol_no!=mol_no)
            {
                (*no_of_molecules)++;
                old_mol_no=mol_no;
            }
            if (strcmp(atom_name,"PL")==0)             //For the P atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].P[i]=coordinate[i];
                //Increasing the count of HPO molecules
                //(*no_of_molecules)++;
            }
            else if (!(OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].OHL_1[i]=coordinate[i];
                OHL_read=1;
            }
            else if (!(HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].HOL_1[i]=coordinate[i];
                HOL_read=1;
            }
            else if ((OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].OHL_2[i]=coordinate[i];
                OHL_read=0;
            }
            else if ((HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].HOL_2[i]=coordinate[i];
                HOL_read=0;
            }
            else if (!(O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].O2L_1[i]=coordinate[i];
                O2L_read=1;
            }
            else if ((O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].O2L_2[i]=coordinate[i];
                O2L_read=0;
            }
        
        }
    }
    (*no_of_molecules)++;
    (*no_of_molecules)=(*no_of_molecules)/(*conf_number);
    fclose(fp_in);
}

void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,char ** atom_name_list, int no_of_molecules)
{
    int i,j,k;
    //Input K molecules
    for(i=0;i<no_of_molecules;i++){
        for(j=0;j<3;j++){
            Kmolecules[i].posn[j]=x[i][j]*10;
        }
    }
    int atom_no=no_of_molecules;
    int index=0;
    int OHL_read=0,O2L_read=0,HOL_read=0;
    //Input HPO molecules
    for(j=0;j<no_of_molecules;j++)
    {
        //printf("molecule %d\n",j);
        for(k=0;k<7;k++)
        {
            if (strcmp(atom_name_list[index%7],"PL")==0)             //For the P atom
            {
                //printf("inside pl\n");
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].P[i]=x[atom_no][i]*10;
                index++;
                atom_no++;
                //Increasing the count of HPO molecules
                //(*no_of_molecules)++;
            }
            else if (!(OHL_read)&&(strcmp(atom_name_list[index%7],"OHL")==0))//For the OHL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].OHL_1[i]=x[atom_no][i]*10;
                OHL_read=1;
                index++;
                atom_no++;
            }
            else if (!(HOL_read)&&(strcmp(atom_name_list[index%7],"HOL")==0))//For the HOL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].HOL_1[i]=x[atom_no][i]*10;
                HOL_read=1;
                index++;
                atom_no++;
            }
            else if ((OHL_read)&&(strcmp(atom_name_list[index%7],"OHL")==0))//For the OHL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].OHL_2[i]=x[atom_no][i]*10;
                OHL_read=0;
                index++;
                atom_no++;
            }
            else if ((HOL_read)&&(strcmp(atom_name_list[index%7],"HOL")==0))//For the HOL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].HOL_2[i]=x[atom_no][i]*10;
                HOL_read=0;
                index++;
                atom_no++;
            }
            else if (!(O2L_read)&&(strcmp(atom_name_list[index%7],"O2L")==0))//For the O2L_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].O2L_1[i]=x[atom_no][i]*10;
                O2L_read=1;
                index++;
                atom_no++;
            }
            else if ((O2L_read)&&(strcmp(atom_name_list[index%7],"O2L")==0))//For the O2L_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[j].O2L_2[i]=x[atom_no][i]*10;
                O2L_read=0;
                index++;
                atom_no++;
            }
        }

    }
}

void XTC_reader(struct t_fileio* fio,FILE* fp_top,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number,real time_to_start)
{
        char **atom_name_list;
        atom_name_list=(char **)malloc(sizeof(char *)*7);
        TOP_reader(fp_top,no_of_molecules,"HPO",atom_name_list);
        int natoms,i,j;
        int64_t step;
        real time,prec;
        matrix box;
        gmx_bool bOK;
        rvec* x;
        int a = read_first_xtc(fio,&natoms,&step,&time,box,&x,&prec,&bOK);
        for(i=0;i<3;i++)
        {
            boxlength[i]=box[i][i]*10;
        }
        *conf_number=0;
        *start_mol_no=*no_of_molecules;
        do
        {
            if(time>=time_to_start)
            {
                molecule_entry(&molecules[(*no_of_molecules)*(*conf_number)],&Kmolecules[(*no_of_molecules)*(*conf_number)],x,atom_name_list,*no_of_molecules);
                (*conf_number)++;
            }

        }
        while(read_next_xtc(fio,natoms,&step,&time,box,x,&prec,&bOK));
        close_xtc(fio);
}


void skip_comments(FILE *fp_top,char *line)
{
    char semi_colon=' ';
    while(1)
    {
        if(fgets(line, LLEN, fp_top) == NULL)
            return;
        sscanf(line,"%c",&semi_colon);
        if(semi_colon!=';')
            break;
    }
}
// void skip_blank_lines(FILE *fp_top,char *line)
// {
//     while(1)
//     {
//         line[0]='\0';
//         if(fgets(line, LLEN, fp_top) == NULL)
//             return;
//         if(line[0]!='\0')   
//             break;
//         else
//             printf("Empty line \n");
//     }
// }

void read_mol_order(FILE *fp_top,char *line,char** atom_list)
{
    char atom_name[LLEN],mol_name[LLEN];
    int num=0,old_num=0;
    int index=0;
    do
    {
        old_num=num;
        sscanf(line,"%d %s %*d %s %*s %*d %*lf",&num,atom_name,mol_name);
        if(old_num==num)
            break;
        else
        {
            //printf("%s\n",atom_name);
            atom_list[index]=(char *)malloc(sizeof(char)*strlen(atom_name));
            strcpy(atom_list[index++],atom_name);
        } 
        if(fgets(line, LLEN, fp_top) == NULL)
            return;
    }
    while(1);
}


void TOP_reader(FILE* fp_top,int *no_of_molecules,char *molecule_to_search,char ** atom_name_list)
{
    char line[LLEN],atom_name[LLEN],mol_name[LLEN],tag[LLEN];
    int molecule_found=0;
    while(1)
    {
        tag[0]='\0';
        if(fgets(line, LLEN, fp_top) == NULL)
            return;
        sscanf(line,"[ %s ]",tag);
        if(strcasecmp(tag,"moleculetype")==0)
        {
            //Skips comments which start with ;
            skip_comments(fp_top,line);
            sscanf(line,"%s",mol_name);
            if(strcasecmp(molecule_to_search,mol_name)==0)
            {
                molecule_found=1;
            }
        }
        else if((strcasecmp(tag,"atoms")==0)&&(molecule_found==1))
        {
            skip_comments(fp_top,line);
            read_mol_order(fp_top,line,atom_name_list);
        }
        else if((strcasecmp(tag,"molecules")==0))
        {
            skip_comments(fp_top,line);
            do
            {
                sscanf(line,"%s %d",mol_name,no_of_molecules);
                if(strcasecmp(molecule_to_search,mol_name)==0)
                {
                    //printf("Number of molecules: %d\n",*no_of_molecules);
                    break;
                }
                if(fgets(line, LLEN, fp_top) == NULL)
                    return;   
            } while (1);
        }
        //printf("%s",tag);
    }
    fclose(fp_top);
}
