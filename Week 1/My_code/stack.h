#include <stdio.h>
//#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>

typedef struct node
{
    int data;
    struct node* next;
}node;

typedef struct
{
    node* top;
}stack;

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


