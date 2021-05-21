#include <stdio.h>
//#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

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
}stack;

//Creates a node and adds it to stack. Its value is provided as an argument
void add_node_given_value(stack existing_stack, int val)
{
    node* p=malloc(sizeof(node));
    p->data=val;
    if(existing_stack.top==NULL)
    {
        existing_stack.top=p;
        p->next=NULL;
    }
    else
    {
        p->next=existing_stack.top;
        existing_stack.top=p;
    }
}

//Adds a node to the stack. The pointer to the node is provided as an argument
void add_node(stack existing_stack, node* p)
{
    if(existing_stack.top==NULL)
    {
        existing_stack.top=p;
        p->next=NULL;
    }
    else
    {
        p->next=existing_stack.top;
        existing_stack.top=p;
    }
}

//Pops node from the top of stack, frees the node and returns the value of interger
int pop_and_return_value(stack existing_stack)
{
    node *temp;
    if(existing_stack.top==NULL)
        return -1;
    temp=existing_stack.top;
    existing_stack.top=existing_stack.top->next;
    //Consider making c a global variable to avoid repeated declaration
    int c=temp->data;
    free(temp);
    return c;
}

//Pops elements from stack A and pushes them onto stack B until stack A is empty
void pop_from_a_into_b(stack a,stack b)
{
    node *temp;
    while(a.top!=NULL)
    {
        temp=a.top;
        a.top=a.top->next;
        add_node(b,temp);
    }
}


