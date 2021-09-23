#include "data_manip.h"

vector<int> list_to_one_hot(vector<int> list , int max_num)
{
    vector<int> one_hot(0,max_num);
    for(auto x: list)
        one_hot[x]=1;
    return one_hot;
}
vector<int> list_to_one_hot(stack list, int max_num)
{
    vector<int> one_hot(0,max_num);
    node* temp=list.top;
    while (temp!=NULL)
    {
        one_hot[temp->data]=1;
        temp=temp->next;
    }
    return one_hot;
}
vector<int> one_hot_to_list(vector<int> one_hot)
{
    vector<int> list;
    for(int i=0;i<one_hot.size();i++)
        if(one_hot[i]==1)
            list.push_back(i);
    return list;
}