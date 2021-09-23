#ifndef DATAMANIP_H
#define DATAMANIP_H

#include <vector>
//#include <stack>
#include "stack.h"

using namespace std;

//#define LLEN 300
vector<int> list_to_one_hot(vector<int> list , int max_num);
vector<int> list_to_one_hot(stack list, int max_num);
vector<int> one_hot_to_list(vector<int> one_hot);

//int parallelism_enabled=1;

#endif