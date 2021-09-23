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
                    molecules[*no_of_molecules].posn[PL][i]=coordinate[i];
                //Increasing the count of HPO molecules
                //(*no_of_molecules)++;
            }
            else if (!(OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[OHL_1][i]=coordinate[i];
                OHL_read=1;
            }
            else if (!(HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[HOL_1][i]=coordinate[i];
                HOL_read=1;
            }
            else if ((OHL_read)&&(strcmp(atom_name,"OHL")==0))//For the OHL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[OHL_2][i]=coordinate[i];
                OHL_read=0;
            }
            else if ((HOL_read)&&(strcmp(atom_name,"HOL")==0))//For the HOL_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[HOL_2][i]=coordinate[i];
                HOL_read=0;
            }
            else if (!(O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_1 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[O2L_1][i]=coordinate[i];
                O2L_read=1;
            }
            else if ((O2L_read)&&(strcmp(atom_name,"O2L")==0))//For the O2L_2 atom
            {
                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                for(i=0;i<3;i++)
                    molecules[*no_of_molecules].posn[O2L_2][i]=coordinate[i];
                O2L_read=0;
            }
        
        }
    }
    (*no_of_molecules)++;
    (*no_of_molecules)=(*no_of_molecules)/(*conf_number);
    fclose(fp_in);
}


void molecule_entry(HPO molecules[],K Kmolecules[],rvec* x,int * atom_index, int no_of_molecules)
{
    int i,j,k;
    //Input K molecules
    for(i=0;i<no_of_molecules;i++){
        for(j=0;j<DIM;j++){
            Kmolecules[i].posn[j]=x[i][j]*10;
        }
    }
    int atom_no=no_of_molecules;
    int index=0;
    //Input HPO molecules
    for(j=0;j<no_of_molecules;j++)
    {
        //printf("molecule %d\n",j);
        for(k=0;k<HPO_ATOM_COUNT;k++)
        {

            //printf("inside pl\n");
            //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
            for(i=0;i<DIM;i++)
                molecules[j].posn[atom_index[index%HPO_ATOM_COUNT]][i]=x[atom_no][i]*10;
            index++;
            atom_no++;
            //Increasing the count of HPO molecules
            //(*no_of_molecules)++;
        }

    }
}

void XTC_reader(struct t_fileio* fio,FILE* fp_top,HPO molecules[],K Kmolecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no,int *conf_number,real time_to_start)
{
        char **atom_name_list;
        atom_name_list=(char **)malloc(sizeof(char *)*HPO_ATOM_COUNT);
        orderOfMolsArray order;
        int molSumCount=0;
        long long count=0;
        long long Kcount=0;
        long long xcount=0;
        moleculeInfoList molInfo;
        TOP_parser(fp_top,&molInfo,&order);
        //printf("%ld %ld\n",molInfo.size(),order.size());
        *no_of_molecules=0;
        for(auto x:order)
        {
            //printf("hi\n");
            if(strcmp(x.name.c_str(),"HPO")==0)
                *no_of_molecules+=x.quantity;
        }
        //printf("no of mols=%d\n",*no_of_molecules);
        //TOP_reader(fp_top,no_of_molecules,"HPO",atom_name_list, &order);
        // for(auto x: order)
        // {
        //     printf("%s %d\n",x.name.c_str(),x.quantity);
        // }
        vector<int> atom_index;
        //int atom_index[HPO_ATOM_COUNT];
        for(auto x: molInfo)
            if(strcmp(x.moleculeName.c_str(),"HPO")==0)
                atom_name_to_index(x.atomName, &atom_index);
        // for(auto x: atom_index)
        //     printf("%d \n",x);
        int natoms,i,j,k;
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
        //printf("hi\n");
        do
        {
            if(time>=time_to_start)
            {
                //int old_size=0;
                //Cofiguration starts
                xcount=0;
                for(auto o: order)
                {
                    if(strcmp(o.name.c_str(),"K")==0)
                    {
                        //printf("hiK\n");
                        //old_size=o.quantity;
                        //continue;
                        for(i=0;i<o.quantity;i++){
                            for(j=0;j<DIM;j++){
                                Kmolecules[i+Kcount].posn[j]=x[i+xcount][j]*10;
                            }
                        }
                        Kcount+=o.quantity;
                        xcount+=o.quantity;
                    }
                    else if(strcmp(o.name.c_str(),"HPO")==0)
                    {
                        //printf("hiHPO\n");
                        //printf("reached here somehow\n");
                        //molecule_entry(&molecules[count],&Kmolecules[Kcount],&(x[xcount]),atom_index,o.quantity);
                        //printf("%d\n",o.quantity);
                        
                        for(j=0;j<o.quantity;j++)
                        {
                            //printf("molecule %d\n",j);
                            int index=0;
                            for(k=0;k<HPO_ATOM_COUNT;k++)
                            {

                                //printf("inside pl\n");
                                //printf("%d %s %s %d %lf %lf %lf\n",atom_no,atom_name,mol_name,mol_no,coordinate[0],coordinate[1],coordinate[2]);
                                
                                for(i=0;i<DIM;i++)
                                    molecules[count+j].posn[atom_index[index]][i]=x[xcount+j*HPO_ATOM_COUNT+k][i]*10;
                                index++;
                                //index%=HPO_ATOM_COUNT;
                                //atom_no++;
                                //Increasing the count of HPO molecules
                                //(*no_of_molecules)++;
                            }

                        }
                        //molSumCount+=o.quantity;
                        count+=o.quantity;
                        //Kcount+=o.quantity;
                        xcount+=o.quantity*(HPO_ATOM_COUNT);
                    }
                    //else break;
                }

                (*conf_number)++;
                //xcount=0;
                //*no_of_molecules=molSumCount;
                //molSumCount=0;

            }

        }

        while(read_next_xtc(fio,natoms,&step,&time,box,x,&prec,&bOK));
        *start_mol_no=*no_of_molecules;
        // for(int i=0;i<*no_of_molecules;i++)
        //     print_HPO(&(molecules[i]));
        // for(int i=0;i<*no_of_molecules;i++)
        //    printf("%lf %lf %lf\n",Kmolecules[i].posn[0],Kmolecules[i].posn[1],Kmolecules[i].posn[2]);
        //printf("%ld %ld %ld\n",count, Kcount, xcount);
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
        sscanf(line,"%d %*s %*d %s %s %*d %*lf",&num,mol_name,atom_name);
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

void TOP_parser(FILE* fp_top, moleculeInfoList* molInfo ,orderOfMolsArray* order)
{
    char line[LLEN],atom_name[LLEN],mol_name[LLEN],tag[LLEN];
    int enterMoleculeType=0;
    int enterAtoms=0;
    int enterMolecules=0;
    int molNum=-1;
    //int molinfoIndex=-1;
    while(1)
    {
        tag[0]='\0';
        if(fgets(line, LLEN, fp_top) == NULL)
            break;
        else if((strcmp(line,"\n")==0)||(strcmp(line,"\r\n")==0))
        {
            enterMolecules=0;
            enterAtoms=0;
            enterMolecules=0;
            continue;
        }
        else if(line[0]==';')
        {
            //printf("encountered comment\n");
            continue;
        } //Encountered comment
        else if(enterMoleculeType==1)
        {
            //printf("entered Moleculetype\n");
            sscanf(line,"%s %*d",mol_name);
            moleculeInfo temp;
            temp.moleculeName=mol_name;
            molInfo->push_back(temp);
            enterMoleculeType=0;
            continue;
        }
        else if(enterAtoms==1)
        {
            //printf("entered atoms\n");
            // if((strcmp(line,"\n")==0)||(strcmp(line,"\r\n")==0))
            // {
            //     printf("end of atoms\n");
            //     enterAtoms=0;
            //     continue;
            // }
            sscanf(line,"%*d %*s %*d %s %s %*d %*lf",mol_name,atom_name);
            int index=0;
            for(auto x: *molInfo)
            {
                if(strcmp(x.moleculeName.c_str(),mol_name)==0)
                {
                    //printf("found\n");
                    break;
                }
                index++;
            }
            if(index<molInfo->size())
                (*molInfo)[index].atomName.push_back(atom_name);

            //printf("size: %ld\n",molInfo[index].atomName.size());
            continue;

        }
        else if(enterMolecules==1)
        {
            //printf("entered molecules\n");
            // if(line[0]=='\0')
            // {
            //     enterMolecules=0;
            //     continue;
            // }
            //printf("hi\n");
            sscanf(line,"%s %d",mol_name,&molNum);
            //printf("hi\n");
            orderOfMols temp;
            temp.name=mol_name;
            temp.quantity=molNum;
            //printf("hi\n");
            order->push_back(temp); 
            //printf("Order size: %ld\n", order.size());

        }

        sscanf(line,"[ %s ]",tag);

        if(strcasecmp(tag,"moleculetype")==0)
        {
            enterMoleculeType=1;
            continue;
        }
        else if(strcasecmp(tag,"atoms")==0)
        {
            enterAtoms=1;
            continue;
        }
        else if(strcasecmp(tag,"molecules")==0)
        {
            enterMolecules=1;
            continue;
        }
        
        //printf("%s",tag);
    }
    //printf("Order size: %ld\n", order.size());
    //printf("MolInfo size: %ld\n", molInfo.size());
    fclose(fp_top);
    // for(auto x:*order)
    // {
    //     printf("%s %d\n",x.name.c_str(),x.quantity);
    // }
    // for(auto x:*molInfo)
    // {
    //     printf("%s - %ld -", x.moleculeName.c_str(),x.atomName.size());
    //     for(auto y:x.atomName)
    //     {
    //         printf("%s ", y.c_str());
    //         printf(" 1 ");
    //     }
    //     printf("\n");
    // }
}


void TOP_reader(FILE* fp_top,int *no_of_molecules,char *molecule_to_search,char ** atom_name_list, orderOfMolsArray *order)
{
    char line[LLEN],atom_name[LLEN],mol_name[LLEN],tag[LLEN];
    int molecule_found=0;
    int actualno=0;
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
                orderOfMols temp;
                temp.name=mol_name;
                temp.quantity=*no_of_molecules;
                order->push_back(temp);

                // if(strcasecmp(molecule_to_search,mol_name)==0)
                // {
                //     //printf("Number of molecules: %d\n",*no_of_molecules);
                //     break;
                // }
                if(fgets(line, LLEN, fp_top) == NULL)
                    break;
                    //return;   
            } while (1);
            //printf("Order size: %ld\n", order->size());
        }
        
        //printf("%s",tag);
    }
    
    fclose(fp_top);
}

void atom_name_to_index(char ** atom_name_list, int* atom_index)
{
    int i;
    int OHL_read=0,O2L_read=0,HOL_read=0;
    for(i=0;i<HPO_ATOM_COUNT;i++)
    {
        if(strcmp(atom_name_list[i],"PL")==0)
            atom_index[i]=PL;
        else if(strcmp(atom_name_list[i],"OHL")==0)
        {
            if(!OHL_read)
            {
                OHL_read=1;
                atom_index[i]=OHL_1;
            }
            else
                atom_index[i]=OHL_2;
        }
        else if(strcmp(atom_name_list[i],"HOL")==0)
        {
            if(!HOL_read)
            {
                HOL_read=1;
                atom_index[i]=HOL_1;
            }
            else
                atom_index[i]=HOL_2;
        }
        else if(strcmp(atom_name_list[i],"O2L")==0)
        {
            if(!O2L_read)
            {
                O2L_read=1;
                atom_index[i]=O2L_1;
            }
            else
                atom_index[i]=O2L_2;
        }
    }
}

void atom_name_to_index(vector<string> atom_name_list, vector<int> *atom_index)
{
    int i;
    int OHL_read=0,O2L_read=0,HOL_read=0;
    atom_index->resize(HPO_ATOM_COUNT);
    for(i=0;i<HPO_ATOM_COUNT;i++)
    {
        if(strcmp(atom_name_list[i].c_str(),"PL")==0)
            (*atom_index)[i]=PL;
        else if(strcmp(atom_name_list[i].c_str(),"OHL")==0)
        {
            if(!OHL_read)
            {
                OHL_read=1;
                (*atom_index)[i]=OHL_1;
            }
            else
                (*atom_index)[i]=OHL_2;
        }
        else if(strcmp(atom_name_list[i].c_str(),"HOL")==0)
        {
            if(!HOL_read)
            {
                HOL_read=1;
                (*atom_index)[i]=HOL_1;
            }
            else
                (*atom_index)[i]=HOL_2;
        }
        else if(strcmp(atom_name_list[i].c_str(),"O2L")==0)
        {
            if(!O2L_read)
            {
                O2L_read=1;
                (*atom_index)[i]=O2L_1;
            }
            else
                (*atom_index)[i]=O2L_2;
        }
    }
}
