#include "input.h"
#include "stack.h"

//Populates the array of Molecules
void populator(FILE* fp_in,HPO molecules[],coordinates boxlength,int *no_of_molecules,int *start_mol_no)
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
    printf("BoxLength = %lf %lf %lf \n",boxlength[0],boxlength[1],boxlength[2]);
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
}

//Checks and reports whether the strict definition of connectedness results in the same or lesser connections
void strict_vs_relaxed(HPO molecules[],coordinates boxlength,int no_of_molecules)
{
    int i,j;
    int strict=0,not_strict=0;
    for(i=0;i<no_of_molecules-1;i++)
        for(j=i+1;j<no_of_molecules;j++)
        {
            //printf("Are molecules %d and %d connected: %d\n",i,j,connected_molecules(molecules[i],molecules[j],boxlength));
            strict+=connected_molecules_strict(molecules[i],molecules[j],boxlength);
            not_strict+=connected_molecules(molecules[i],molecules[j],boxlength);
        }
    printf("Strict : %d  Not Strict: %d\n",strict,not_strict);
}


//Constructs adjacency list from the array of molecules by checking connectednes between molecules
void adjacency_list_constructor(HPO molecules[],coordinates boxlength,int no_of_molecules,stack adjacency_list[])
{
    int i,j;
    for(i=0;i<no_of_molecules;i++)
        adjacency_list[i].top=NULL;
    for(i=0;i<no_of_molecules-1;i++)
    {
        for(j=i+1;j<no_of_molecules;j++)
        {
            if(connected_molecules(molecules[i],molecules[j],boxlength))
            {
                add_node_given_value(&adjacency_list[i],j);
                add_node_given_value(&adjacency_list[j],i);

                //Shows which molecules are connected
                printf("Nodes %d and %d are connected\n",i,j);
                //print_stack(&adjacency_list[i]);
                //print_stack(&adjacency_list[j]);
            }
        }
    }
}

void dfs_util(stack* to_search,stack adjacency_list[],int visited[],int no_of_molecules,int cluster_number)
{
    int c;
    while(to_search->top!=NULL)
    {
        c=pop_and_return_value(to_search);
        visited[c]=cluster_number;
        pop_after_checking_visited(&adjacency_list[c],to_search,visited);
    }
}

void dfs(stack adjacency_list[],int visited[],int no_of_molecules)
{
    stack to_search;
    to_search.top=NULL;
    int i=0;
    int cluster_number=0;
    while(i<no_of_molecules)
    {
        if(visited[i]==-1)
        {
            visited[i]=cluster_number;
            if(adjacency_list[i].top!=NULL)
            {
                //print_stack(&adjacency_list[i]);
                //add_node_given_value(&to_search,i);
                //print_stack(&to_search);
                pop_after_checking_visited(&adjacency_list[i],&to_search,visited);
                //printf("Popped after checking visited\n");
                //print_stack(&adjacency_list[i]);
                //print_stack(&to_search);
                dfs_util(&to_search, adjacency_list,visited,no_of_molecules,cluster_number);
            }
            cluster_number++;
        }
        i++;
    }
}



int main(int argc,char *argv[])
{
    FILE *fp_in, *fp_out;
    char filename[NAME], xname[2][NAME], line[LLEN];
    char atom_name[LLEN],mol_name[LLEN];
    int atom_no,mol_no;
    int start_mol_no=-1;
    int no_of_molecules=0;
    int i;
    double box[3][2]={0.0}, temp;
    static HPO molecules[MAX_M]; //only for 5000 molecules
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
                {
                    printf("cannot open Infile \n"); 
                    exit(0); 
                }
                break;
            case 'o':
                if(( fp_out  =fopen(optarg,"w"))==NULL)
                {
                    printf("cannot open Outfile \n"); 
                    exit(0); 
                }
                break;
        }
    }
    /*------------------------------ read the arguments-----------------------------------------*/

    //Populate the array
    populator(fp_in,molecules,boxlength,&no_of_molecules,&start_mol_no);

    //Check if strictness matters
    strict_vs_relaxed(molecules,boxlength,no_of_molecules);

    //int end_mol_no=start_mol_no+no_of_molecules;


    //Prints all the molecules
    // printf("%d %d",start_mol_no,end_mol_no);
    // for(i=start_mol_no;i<end_mol_no;i++)
    //     print_HPO(molecules[i]);

    //int connectedness[MAX_M]={0};

    stack adjacency_list[no_of_molecules];
   // int visited[no_of_molecules]={-1};
    adjacency_list_constructor(molecules,boxlength,no_of_molecules,adjacency_list);
    
    //Printing the adjacency list
    for(i=0;i<no_of_molecules;i++)
    {
        printf("%d :",i);
        print_stack(&adjacency_list[i]);
    }


    int visited[no_of_molecules];
    for(i=0;i<no_of_molecules;i++)
        visited[i]=-1;


    dfs(adjacency_list,visited,no_of_molecules);

    //Printing the cluster number
    for(i=0;i<no_of_molecules;i++)
        printf("%d molecule belongs to cluster number %d\n",i+start_mol_no,visited[i]);

    /*If the adjacency lists are empty after the operation, it is an indication that
    all the vertices were processed during DFS*/    
    for(i=0;i<no_of_molecules;i++)
    {
        printf("%d :",i);
        print_stack(&adjacency_list[i]);
    }

    return 0;
}