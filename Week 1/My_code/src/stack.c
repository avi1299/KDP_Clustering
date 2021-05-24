#include "stack.h"

void print_node(node a)
{
    printf("%d  %p\n",a.data,a.next);
}

void print_stack(stack *a)
{
    node *temp;
    temp=a->top;
    while(temp!=NULL)
    {
        printf("%d <- ",temp->data);
        temp=temp->next;
    }
    printf("NULL\n");
}

void print_stack_elements(stack* a)
{
    node *temp;
    temp=a->top;
    while(temp!=NULL)
    {
        printf("%d ",temp->data);
        temp=temp->next;
    }
    printf("\n");
}

//Creates a node and adds it to stack. Its value is provided as an argument
void add_node_given_value(stack *existing_stack, int val)
{
    node* p=(node*)malloc(sizeof(node));
    p->data=val;
    if(existing_stack->top==NULL)
        p->next=NULL;
    else
        p->next=existing_stack->top;
    existing_stack->top=p;
    existing_stack->length++;

}

//Adds a node to the stack. The pointer to the node is provided as an argument
void add_node(stack *existing_stack, node* p)
{
    if(existing_stack->top==NULL)
        p->next=NULL;
    else
        p->next=existing_stack->top;
    existing_stack->top=p;
    existing_stack->length++;

}

//Pops node from the top of stack, frees the node and returns the value of interger
int pop_and_return_value(stack *existing_stack)
{
    node *temp;
    if(existing_stack->top==NULL)
        return -1;
    temp=existing_stack->top;
    existing_stack->top=(existing_stack->top)->next;
    existing_stack->length--;
    //Consider making c a global variable to avoid repeated declaration
    int c=temp->data;
    free(temp);
    return c;
}

//Pops node from top of stack and returns pointer to node
node* pop_and_return_node(stack *existing_stack)
{
    node *temp;
    if(existing_stack->top==NULL)
        return NULL;
    temp=existing_stack->top;
    existing_stack->top=(existing_stack->top)->next;
    existing_stack->length--;
    return temp;
}

//Pops elements from stack A and pushes them onto stack B until stack A is empty
void pop_from_a_into_b(stack *a,stack *b)
{
    node *temp;
    while(a->top!=NULL)
    {
        temp=a->top;
        a->top=(a->top)->next;
        a->length--;
        add_node(b,temp);
    }
}

//Pops elements from stack A and pushes them onto stack B until stack A is empty. Checks if the elements have been visited or not
void pop_after_checking_visited(stack *a,stack *b,int visited[])
{
    node *temp;
    while(a->top!=NULL)
    {
        temp=a->top;
        a->top=(a->top)->next;
        a->length--;
        if(visited[temp->data]==-1)
            add_node(b,temp);
    }
}

