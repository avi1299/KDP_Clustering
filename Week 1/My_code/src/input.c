#include "input.h"


//Assumptions: All the molecules are listed together and in order of molecules P OHL HOL OHL HOL O2L O2L
void PDB_reader(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no)
{
    int i,atom_no,mol_no;
    char line[LLEN],atom_name[LLEN],mol_name[LLEN];
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
        sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
        //Executes below when it encounters HPO molecules. It extracts information for every atom in every HPO molecule 
        if(strcmp(mol_name,"HPO")==0)
        {
            //Increasing the count of HPO molecules
            (*no_of_molecules)++;
            //Noting where the first HPO molecule was encountered
            if(*start_mol_no==-1)
                *start_mol_no=mol_no;
            //For the P atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].P[i]=coordinate[i];
            }
            //For the OHL_1 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_1[i]=coordinate[i];
            }
            //For the HOL_1 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_1[i]=coordinate[i];
            }
            //For the OHL_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_2[i]=coordinate[i];
            }
            //For the HOL_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_2[i]=coordinate[i];
            }
            //For the O2L_1 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_1[i]=coordinate[i];
            }
            //For the O2L_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_2[i]=coordinate[i];
            }
        
        }
    }
    fclose(fp_in);
}

//Assumptions: All the molecules are listed together and in order of molecules OHL HOL P O2L OHL HOL O2L
void PDB_reader_alternate(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no)
{
    int i,atom_no,mol_no;
    char line[LLEN],atom_name[LLEN],mol_name[LLEN];
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
        sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
        //Executes below when it encounters HPO molecules. It extracts information for every atom in every HPO molecule 
        if(strcmp(mol_name,"HPO")==0)
        {
            //Increasing the count of HPO molecules
            (*no_of_molecules)++;
            //Noting where the first HPO molecule was encountered
            if(*start_mol_no==-1)
                *start_mol_no=mol_no;
            //For the OHL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_1[i]=coordinate[i];

            }
            //For the HOL_1 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_1[i]=coordinate[i];
            }            
            //For the P atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].P[i]=coordinate[i];
            }
            //For the O2L_1 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_1[i]=coordinate[i];
            }
            //For the OHL_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_2[i]=coordinate[i];
            }
            //For the HOL_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_2[i]=coordinate[i];
            }
            //For the O2L_2 atom
            {
                fgets(line, LLEN, fp_in);
                sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_2[i]=coordinate[i];
            }
        
        }
    }
    fclose(fp_in);
}

//Assumption all the molecules are listed together and though order of molecules doesn't matter, all the 7 atoms corresponding to a molecule should be together
//i.e. one after the other.
void PDB_reader_universal(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no)
{
    int i,atom_no,mol_no;
    char line[LLEN],atom_name[LLEN],mol_name[LLEN];
    int HOL_read=0;
    int OHL_read=0;
    int O2L_read=0;
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
        sscanf(line,"%*s %d %s %s X %d %lf %lf %lf %*lf %*lf",&atom_no,atom_name,mol_name,&mol_no,&coordinate[0],&coordinate[1],&coordinate[2]);
        //Executes below when it encounters HPO molecules. It extracts information for every atom in every HPO molecule 
        if(strcmp(mol_name,"HPO")==0)
        {
            //Noting where the first HPO molecule was encountered
            if(*start_mol_no==-1)
                *start_mol_no=mol_no;
            if (strcmp(atom_name,"PL")==0)             //For the P atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].P[i]=coordinate[i];
                //Increasing the count of HPO molecules
                (*no_of_molecules)++;
            }
            else if (!(OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_1[i]=coordinate[i];
                OHL_read=1;
            }
            else if (!(HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_1[i]=coordinate[i];
                HOL_read=1;
            }
            else if ((OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].OHL_2[i]=coordinate[i];
                OHL_read=0;
            }
            else if ((HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].HOL_2[i]=coordinate[i];
                HOL_read=0;
            }
            else if (!(O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_1[i]=coordinate[i];
                O2L_read=1;
            }
            else if ((O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[mol_no-*start_mol_no].O2L_2[i]=coordinate[i];
                O2L_read=0;
            }
        
        }
    }
    fclose(fp_in);
}