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

void print_node(node a);

void print_stack(stack *a);

void print_stack_elements(stack* a);

//Creates a node and adds it to stack. Its value is provided as an argument
void add_node_given_value(stack *existing_stack, int val);

//Adds a node to the stack. The pointer to the node is provided as an argument
void add_node(stack *existing_stack, node* p);

//Pops node from the top of stack, frees the node and returns the value of interger
int pop_and_return_value(stack *existing_stack);

//Pops node from top of stack and returns pointer to node
node* pop_and_return_node(stack *existing_stack);

//Pops elements from stack A and pushes them onto stack B until stack A is empty
void pop_from_a_into_b(stack *a,stack *b);

//Pops elements from stack A and pushes them onto stack B until stack A is empty. Checks if the elements have been visited or not
void pop_after_checking_visited(stack *a,stack *b,int visited[]);

#endif

