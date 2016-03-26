#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "parser.h"
#include "hashtable.h"
#include "mna.h"

int size;
element *tmp2 = NULL;

void init_A_array(){
  
  //initialize two-dimensional A array
  A_array = (double*)calloc((size*size), sizeof(double));
  
}

void init_C_array(){
  
  //initialize two-dimensional A array
  C_array = (double*)calloc((size*size), sizeof(double));
  
}

void init_B_array(){
  
  //initialize one-dimensional B array
  B_array = (double*)calloc(size, sizeof(double));
}

void print_A_array(){
  int i, j;
  
  printf("\nThe %dx%d A array is :\n\n", size, size);
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      printf(" %.4lf \t", A_array[(j*size) + i]);
    }
    printf("\n");
  }
  printf("\n\n");
}


void print_C_array(){
  int i, j;
  
  printf("\nThe %dx%d C array is :\n\n", size, size);
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      printf(" %.4lf \t", C_array[(j*size) + i]);
    }
    printf("\n");
  }
  printf("\n\n");
}


void print_B_array(){
  int i;
  
  printf("\nThe %dx1 B array is:\n\n", size);
  
  for(i=0; i < size; i++){
    printf("%.4lf\n", B_array[i]);
  }
  printf("\n\n");
}

void mna_tran(){
  int pos, neg;
  int k = 0;

  
  init_C_array();
  
  tmp2 = head;
  
  while(tmp2!=NULL){
    if(tmp2->type_id== 'v'){
      k++;
    }
    //take an iductance element
    else if(tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(tmp2->type);
      
      //put to -L piece of C_tilda matrix the -l_k value
      C_array[(n+k)*size +(n+k)] -= temp->value;
      
      k++;
    }
    //take a volt source element
    else if(tmp2->type_id== 'c'){
      
      Capacity* temp = (Capacity*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	C_array[(pos*size) + pos] += temp->value;	//build the C_array
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	C_array[(neg*size) + neg] += temp->value;	//build the C_array
      }
      
      if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	
	C_array[(neg*size) + pos] -= temp->value;	//build the C_array
	C_array[(pos*size) + neg] -= temp->value;	//build the C_array
      
      }
    
    }
    
    tmp2=tmp2->next;
  }
  
}


void mna(int n, int m2){
  int k = 0;
  int pos, neg;

  
  size = n + m2;	//take the size of the A and B array
  
  init_A_array();
  
  init_B_array();
  
  tmp2 = head;
  
  while(tmp2!=NULL){
    
    //take a volt source element
    if(tmp2->type_id== 'v'){
      
      Volt_Source* temp = (Volt_Source*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	A_array[((n + k)*size) + pos] += 1;	//build the A_array
	A_array[(pos*size) + (n + k)] += 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	A_array[((n + k)*size) + neg] -= 1;	//build the A_array
	A_array[(neg*size) + (n + k)] -= 1;
      }
      
      B_array[n+k] += temp->value;	//build the B_array
      
      k++;
    
    }
    //take an iductance element
    else if(tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	A_array[(n + k)*size + pos] += 1;	//build the A_array
	A_array[(pos*size) + (n + k)] += 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	A_array[(n + k)*size + neg] -= 1;	//build the A_array
	A_array[(neg*size) + (n + k)] -= 1;
      }
      k++;
    }
    //take a resistance element
    else if(tmp2->type_id == 'r'){
      
      Resistance* temp = (Resistance*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
      }
      
      //in case that both nodes aren't the ground node
      if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	A_array[(pos*size) + neg] -= (1/(temp->value));		//build the A_array
	A_array[(neg*size) + pos] -= (1/(temp->value));
	A_array[(pos*size) + pos] += (1/(temp->value));
	A_array[(neg*size) + neg] += (1/(temp->value));
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	A_array[(neg*size) + neg] += (1/(temp->value));
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	A_array[(pos*size) + pos] += (1/(temp->value));
      } 
    }
    //take a power source element
    else if(tmp2->type_id == 'i'){
      
      Power_Source* temp = (Power_Source*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	//printf("re\n");
	B_array[pos] -= temp->value;	//build the A_array
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	B_array[neg] += temp->value;	//build the A_array
	//printf("%g\n",temp->value);
      }
      
    }
    
    tmp2=tmp2->next;
  
  }
  
}
