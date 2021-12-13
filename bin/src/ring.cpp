#include "ring.h"


void ringDriver(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], int* ION_list, 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag, int verbose_level, vector<ringCandidate> *CSet, 
    pathArray* CSSSR,vector<ringElements> *CSSSR_Elements)
{
    makePIDmatrix(adjacency_matrix, ION_list, no_of_molecules, D, P, P_dash,(verbose_level>=0));
    ringCandidateSearch(CSet, no_of_molecules, D,P,P_dash);
    //printf("Problem is nnot here\n");
    findSSSR(CSSSR,CSet);

    for(int i=0;i<no_of_molecules;i++)
        for(int j=0;j<no_of_molecules;j++)
        {
            delete P[i][j];
            delete P_dash[i][j];
        }
    //printf("Problem is nnot here\n");
    // for(auto x: CSet)
    //     printf("CNum:%5d | P:%5ld | P_Dash:%5ld\n",x.CNum,x.P->size(),x.P_dash->size());
    ringElements* temp;
    //CSSSR_Elements->clear();
    for(auto x: *CSSSR)
    {
        temp=ringToElements(&x);
        CSSSR_Elements->push_back(*temp);
        //printRingElements(temp);
        delete temp;
    }
    if(verbose_level>=1)
        printRingElementsArray(CSSSR_Elements);
    if(verbose_level>=2)
        printPathArray(CSSSR);

    CSet->clear();
    CSSSR->clear();
}

void makePIDmatrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag)
{
    int i,j,k;
    path s;
    int strong_level=WEAK;
    if(strong_flag)
        strong_level=STRONG;
    #pragma omp parallel for private(s,j)
    for(i=0;i<no_of_molecules;i++)
        for(j=0;j<no_of_molecules;j++)
        {
            P[i][j]= new pathArray;
            P_dash[i][j]= new pathArray;
            if(adjacency_matrix[UNDIRECTED_GRAPH][i][j]>=strong_level)
            {
                D[i][j]=1;
                if(i<j)
                    s.push_back(make_pair(i,j));
                else
                    s.push_back(make_pair(j,i));
                P[i][j]->push_back(s);
                s.clear();
            }
                
            else
                D[i][j]=numeric_limits<int>::max();
            
        }

    // for(i=0;i<no_of_molecules;i++)
    //     for(j=0;j<no_of_molecules;j++)
    //     {
    //         //printf("pathArray from %d to %d:\n", i, j);
    //         printPathArray(P[i][j]);
    //     }

    pathArray* temp;
    for(k=0;k<no_of_molecules;k++)
        for(i=0;i<no_of_molecules;i++)
            for(j=0;j<no_of_molecules;j++)
                if(i!=j&&i!=k&&j!=k&&D[i][k]!=numeric_limits<int>::max()&&D[k][j]!=numeric_limits<int>::max())
                {
                    //printf("hi\n");
                    if(D[i][j]>D[i][k]+D[k][j] && D[i][k]+D[k][j]<=RING_LENGTH_CUTOFF)//The second condition was added to limit the length of the path between any two nodes
                    {
                        //printf("hi1\n");
                        if(D[i][j]==D[i][k]+D[k][j]+1)
                        {
                            *P_dash[i][j]=*P[i][j];
                        }
                        else
                            P_dash[i][j]->clear();
                        D[i][j]=D[i][k]+D[k][j];
                        delete P[i][j];
                        P[i][j]=addPathArray(P[i][k],P[k][j]);
                    }
                    else if(D[i][j]==D[i][k]+D[k][j])
                    {
                        temp=addPathArray(P[i][k],P[k][j]);
                        appendPathArray(P[i][j],temp);
                        delete temp;
                    }
                        
                    else if(D[i][j]==D[i][k]+D[k][j]-1)
                    {
                        temp=addPathArray(P[i][k],P[k][j]);
                        appendPathArray(P_dash[i][j],temp);
                        delete temp;
                    }
                }
    
    // for(i=0;i<no_of_molecules;i++)
    //     for(j=0;j<no_of_molecules;j++)
    //     {
    //         //printf("pathArray from %d to %d:\n", i, j);
    //         printPathArray(P[i][j]);
    //     }


//printf("Made PID\n");
}

void makePIDmatrix(int adjacency_matrix[2][MAX_MOLECULES][MAX_MOLECULES], int* ION_list, 
    int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES],int strong_flag)
{
    int i,j,k;
    path s;
    int strong_level=WEAK;
    if(strong_flag)
        strong_level=STRONG;
    #pragma omp parallel for private(s,j)
    for(i=0;i<no_of_molecules;i++)
        for(j=0;j<no_of_molecules;j++)
        {
            // if(P[i][j]!=NULL)
            // delete P[i][j];
            // if(P_dash[i][j]!=NULL)
            // delete P_dash[i][j];
            P[i][j]= new pathArray;
            P_dash[i][j]= new pathArray;
            if(adjacency_matrix[UNDIRECTED_GRAPH][ION_list[i]][ION_list[j]]>=strong_level)
            {
                D[i][j]=1;
                if(i<j)
                    s.push_back(make_pair(i,j));
                else
                    s.push_back(make_pair(j,i));
                P[i][j]->push_back(s);
                s.clear();
            }
                
            else
                D[i][j]=numeric_limits<int>::max();
            
        }

    // for(i=0;i<no_of_molecules;i++)
    // {
    //     for(j=0;j<no_of_molecules;j++)
    //     {
    //         printf("%d ",D[i][j]);
    //     }
    //     printf("\n");
    // }

    // for(i=0;i<no_of_molecules;i++)
    //     for(j=0;j<no_of_molecules;j++)
    //     {
    //         //printf("pathArray from %d to %d:\n", i, j);
    //         printPathArray(P[i][j]);
    //     }

    pathArray* temp;
    for(k=0;k<no_of_molecules;k++)
        for(i=0;i<no_of_molecules;i++)
            for(j=0;j<no_of_molecules;j++)
                if(i!=j&&i!=k&&j!=k&&D[i][k]!=numeric_limits<int>::max()&&D[k][j]!=numeric_limits<int>::max())
                {
                    //printf("hi\n");
                    if(D[i][j]>D[i][k]+D[k][j])// && D[i][k]+D[k][j]<=RING_LENGTH_CUTOFF)//The second condition was added to limit the length of the path between any two nodes
                    {
                        //printf("hi1\n");
                        if(D[i][j]==D[i][k]+D[k][j]+1)
                        {
                            *P_dash[i][j]=*P[i][j];
                        }
                        else
                            P_dash[i][j]->clear();
                        D[i][j]=D[i][k]+D[k][j];
                        delete P[i][j];
                        P[i][j]=addPathArray(P[i][k],P[k][j]);
                    }
                    else if(D[i][j]==D[i][k]+D[k][j])
                    {
                        temp=addPathArray(P[i][k],P[k][j]);
                        appendPathArray(P[i][j],temp);
                        delete temp;
                    }
                        
                    else if(D[i][j]==D[i][k]+D[k][j]-1)
                    {
                        temp=addPathArray(P[i][k],P[k][j]);
                        appendPathArray(P_dash[i][j],temp);
                        delete temp;
                    }
                }
    
    // for(i=0;i<no_of_molecules;i++)
    //     for(j=0;j<no_of_molecules;j++)
    //     {
    //         //printf("pathArray from %d to %d:\n", i, j);
    //         printPathArray(P[i][j]);
    //     }


//printf("Made PID\n");
}

void ringCandidateSearch(vector<ringCandidate> *CSet, int no_of_molecules, int D[MAX_MOLECULES][MAX_MOLECULES], pathArray *P[MAX_MOLECULES][MAX_MOLECULES],
    pathArray *P_dash[MAX_MOLECULES][MAX_MOLECULES])
{
    int i,j,CNum;
    CSet->clear();
    ringCandidate temp;
    for(i=0;i<no_of_molecules;i++)
        for(j=0;j<no_of_molecules;j++)
            if(i!=j)
            {
                //Continues if no path exists between i and j or if there is only one such path between i and j  
                if (D[i][j] == numeric_limits<int>::max() || ( (P[i][j]->size() == 1) && (P_dash[i][j]->size() == 0)))
                {
                    continue;
                }
                else 
                {
                    //printf("%d\n", D[i][j]);
                    if(P_dash[i][j]->size()!=0)
                        CNum=2*D[i][j]+1;
                    else
                        CNum=2*D[i][j];
                    temp.CNum=CNum;
                    temp.P=P[i][j];
                    temp.P_dash=P_dash[i][j];
                    CSet->push_back(temp);
                }
            }
    // for(auto x: *CSet)
    // {
    //     printf("%d\n", x.CNum);
    //     printPathArray(x.P);
    //     printPathArray(x.P_dash);
    // }


    //printf("Problem is not after deleteing\n");
}

void findSSSR(pathArray* CSSSR, vector<ringCandidate> *CSet)
{
    //TODO: Implement nRingIdx and nSSSR
    
    CSSSR->clear();
    int i, j;
    pathArray *C;
    path temp;
    //printf("CSet size:%d\n",CSet->size());
    for(auto x: *CSet)
    {
        //printf("x.CNum=%d\n",x.CNum);
        if(x.CNum%2==1)
        {
            C=addPathArrayWithoutCommonElements(x.P,x.P_dash);
        }
        else
        {
            // C=new pathArray;
            // for(auto y:*(x.P))
            //     for(auto z:*(x.P))
            //         if(y!=z)
            //         {
            //             temp=y;
            //             temp.insert(temp.end(), z.begin(), z.end());
            //             C->push_back(temp);
            //             temp.clear();
            //         }
            C=addPathArrayWithoutCommonElements(x.P,x.P);
        }
        //printf("should be bad here\n");
        for(auto ring: *C)
            pathArrayXORandAdd(CSSSR, &ring);
        //C->clear();
        delete C;
    }
    

    removeDuplicateRings(CSSSR);

    //printPathArray(CSSSR);

}

void printEdge(edge *intPair)
{
    printf("%d : %d\n", intPair->first, intPair->second);
}
void printPath(path *path)
{
    for(auto x: *path)
        printEdge(&x);

}
void printPathArray(pathArray *patharray)
{
    int count=0;
    for(auto x: *patharray)
    {
        printf("Set %d:\n",count++);
        printPath(&x);
    }
    //printf("\n");
}

void printRingElements(ringElements *ring)
{
    for(auto x: *ring)
        printf("%d ",x);
    printf("\n");
}

void printRingElementsArray(ringElementsArray *ringarray)
{
    int count=0;
    for(auto x: *ringarray)
    {
        printf("Ring %5d : ",count++);
        printRingElements(&x);
    }
}

void appendPathArray(pathArray* arr1, pathArray* arr2)
{
    arr1->insert(arr1->end(),arr2->begin(), arr2->end());
    // for(auto x: *arr1)
    //     removeDuplicateEdgePairs(&x);
}

pathArray *addPathArray(pathArray* arr1, pathArray* arr2)
{
    pathArray *arr=new pathArray;
    path temp;
    for(auto x: *arr1)
    {
        for(auto y: *arr2)
        {
            temp=x;
            temp.insert(temp.end(),y.begin(), y.end());
            arr->push_back(temp);
            temp.clear();
        }
    }
    // for(auto x: *arr)
    //     removeDuplicateEdgePairs(&x);


    return arr;

}



ringElements* ringToElements(path* ring)
{
    ringElements *arr=new ringElements;
    path temp;
    set<int> s;
    for(auto x: *ring)
    {
        s.insert(x.first);
        s.insert(x.second);
    }
    for(auto x : s)
        arr->push_back(x);
    return arr;
}

pathArray *addPathArrayWithoutCommonElements(pathArray* arr1, pathArray* arr2)
{
    pathArray *arr=new pathArray;
    path temp;
    for(auto x: *arr1)
    {
        for(auto y: *arr2)
        {
            if(!pathIntersection(&x,&y))
            {
                temp=x;
                temp.insert(temp.end(),y.begin(), y.end());
                arr->push_back(temp);
                temp.clear();
            }
        }
    }
    // for(auto x: *arr)
    //     removeDuplicateEdgePairs(&x);


    return arr;

}

void pathArrayXORandAdd(pathArray *CSSSR, path *ring)
{
    sort(ring->begin(), ring->end(), comp());
    // printf("Stopped before removing duplicates\n");
    // printpath(ring);
    //removeDuplicateEdgePairs(ring);
    //printpath(ring);
    int insert_flag=1;
    for(auto x: *CSSSR)
    {
        if(x==*ring)
        {
            insert_flag=0;
            break;
        }
    }
    if(insert_flag)
    CSSSR->push_back(*ring);
}


bool pathIntersection(path *small, path *big)
{
    for (auto small_elt: *small) {
        //Keeps finding if elements of small arr in big arr. IF found it returns true else it checks for next value in small until no more values in small.
        if (find(big->begin(), big->end(), small_elt) != big->end()) {
            return true;
        }
    }
    return false;
}

void removeDuplicateRings(pathArray *arr)
{
    int len = arr->size();
    if(len<2)
        return;
    sort(arr->begin(), arr->end(), comppathSizeAsc());
    path small,big;
    for(int i=0;i<arr->size()-1;i++)
    {
        small=(*arr)[i];
        for(int j=i+1;j<arr->size();j++)
        {
            
            big=(*arr)[j];
            if(ringSubset(&small,&big))
            {
                //x->erase(next(x->begin(),i),next(x->begin(),i+2));
                arr->erase(arr->begin()+j);
                //printf("Removed duplicates\n");
                j--;
            }
            
        }
    }
        
}

bool ringSubset(path *small, path *big)
{
    for (auto small_elt: *small) {
        //Keeps finding if elements of small arr in big arr. IF found it checks for next value in small until no more values in small.
        if (find(big->begin(), big->end(), small_elt) == big->end()) {
            return false;
        }
    }
    return true;
}

/*------------------------Older Unnecessary Code--------------------------------*/


void removeDuplicateEdgePairs(path *x)
{
    int len = x->size();
    if(len<2)
        return;
    sort(x->begin(), x->end(), comp());
    
    for(int i=0;i<x->size()-1;i++)
    {
        // edge a,b;
        // a=(*x)[i];
        // b=(*x)[i+1];
        if((*x)[i]==(*x)[i+1])
        {
            //x->erase(next(x->begin(),i),next(x->begin(),i+2));
            x->erase(x->begin()+i);
            x->erase(x->begin()+i);
            //printf("Removed duplicates\n");
            i--;
        }
        
    }
}






