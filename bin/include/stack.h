#ifndef STACK_HEADER
#define STACK_HEADER
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//Structure for node
typedef struct node
{
    int data;
    struct node* next;
}node;

//Structure for stack
typedef struct
{
    node* top;
    int length;
}stack;

/**
 * @brief Prints the node
 * 
 * @param a node
 */
void print_node(node a);

/**
 * @brief Prints the stack
 * 
 * @param a stack*
 */
void print_stack(stack *a);

/**
 * @brief Prints the stack elemets without the NULL or arrows
 * 
 * @param a stack*
 */
void print_stack_elements(stack* a);

/**
 * @brief Creates a node and adds it to stack. Its value is provided as an argument
 * 
 * @param existing_stack stack*
 * @param val int
 */
void add_node_given_value(stack *existing_stack, int val);

/**
 * @brief Adds a node to the stack. The pointer to the node is provided as an argument
 * 
 * @param existing_stack stack*
 * @param p node*
 */
void add_node(stack *existing_stack, node* p);

/**
 * @brief Pops node from the top of stack, frees the node and returns the value of interger
 * 
 * @param existing_stack stack *
 * @return int 
 */
int pop_and_return_value(stack *existing_stack);

/**
 * @brief Pops node from top of stack and returns pointer to node
 * 
 * @param existing_stack stack* 
 * @return node* 
 */
node* pop_and_return_node(stack *existing_stack);

/**
 * @brief Pops elements from stack A and pushes them onto stack B until stack A is empty

 * 
 * @param a stack* 
 * @param b stack* 
 */
void pop_from_a_into_b(stack *a,stack *b);

/**
 * @brief Pops elements from stack A and pushes them onto stack B until stack A is empty.
 * Checks if the elements have been visited or not
 * 
 * @param a stack*
 * @param b stack*
 * @param visited int*
 */
void pop_after_checking_visited(stack *a,stack *b,int visited[]);

/**
 * @brief Pops all elements of the stack and return nothing
 * 
 * @param a stack*
 */
void empty_stack(stack *a);

#endif

