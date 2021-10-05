#ifndef DATAMANIP_H
#define DATAMANIP_H

#include <vector>
//#include <stack>
#include "stack.h"

using namespace std;

//#define LLEN 300

/**
 * @brief Converts a vector of integers to one-hot encoded vector
 * 
 * @param list vector<int>
 * @param max_num int
 * @return vector<int> 
 */
vector<int> list_to_one_hot(vector<int> list , int max_num);

/**
 * @brief Converts a stack of integers to one-hot encoded vector
 * 
 * @param list stack
 * @param max_num int
 * @return vector<int> 
 */
vector<int> list_to_one_hot(stack list, int max_num);

/**
 * @brief Converts a one-hot encoded vector to a vector of integers
 * 
 * @param one_hot vector<int>
 * @return vector<int> 
 */
vector<int> one_hot_to_list(vector<int> one_hot);


#endif