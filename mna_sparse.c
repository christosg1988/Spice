#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "parser.h"
#include "hashtable.h"
#include "mna_sparse.h"
#include "mna.h"
#include <string.h>


int size;

void mna_sparse(int n, int m2, int n_z){
  int pos, neg, k=0, l=-1;
  
  
  size = n + m2;

  //initialize one-dimensional B array
  B_array_cs = (double*)calloc(size, sizeof(double));
  
  Asparse = cs_spalloc(size, size, n_z , 1, 1);
  Asparse->nz = n_z;

  element *tmp2 = NULL;

  tmp2 = head;

  while(tmp2!=NULL){

    //take a volt source element
    if(tmp2->type_id== 'v'){
      
      Volt_Source* temp = (Volt_Source*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	
	//build triplet
	l++;
	
	Asparse->i[l] = n + k;
	Asparse->p[l] = pos;
	Asparse->x[l] = 1;
	
	l++;
	
	Asparse->i[l] = pos;
	Asparse->p[l] = n + k;
	Asparse->x[l] = 1;
	
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	//build triplet
	l++;
	
	Asparse->i[l] = n + k;
	Asparse->p[l] = neg;
	Asparse->x[l] = -1;
	
	l++;
	
	Asparse->i[l] = neg;
	Asparse->p[l] = n + k;
	Asparse->x[l] = -1;
	
      }

      B_array_cs[n+k] += temp->value;	//build the B_array
      
      k++;
    }
    //take an iductance element
    else if(tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	//build triplet
	l++;
	
	Asparse->i[l] = n + k;
	Asparse->p[l] = pos;
	Asparse->x[l] = 1;
	
	l++;
	
	Asparse->i[l] = pos;
	Asparse->p[l] = n + k;
	Asparse->x[l] = 1;
	
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	//build triplet
	l++;
	
	Asparse->i[l] = n + k;
	Asparse->p[l] = neg;
	Asparse->x[l] = -1;
	
	l++;
	
	Asparse->i[l] = neg;
	Asparse->p[l] = n + k;
	Asparse->x[l] = -1;
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
	
	//build triplet
	l++;
	
	Asparse->i[l] = pos;
	Asparse->p[l] = neg;
	Asparse->x[l] = -1/temp->value;
	
	l++;
	
	Asparse->i[l] = neg;
	Asparse->p[l] = pos;
	Asparse->x[l] = -1/temp->value;
	
	l++;
	
	Asparse->i[l] = pos;
	Asparse->p[l] = pos;
	Asparse->x[l] = 1/temp->value;
	
	l++;
	
	Asparse->i[l] = neg;
	Asparse->p[l] = neg;
	Asparse->x[l] = 1/temp->value;
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	
	//build triplet
	l++;
	
	Asparse->i[l] = neg;
	Asparse->p[l] = neg;
	Asparse->x[l] = 1/temp->value;
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	
	//build triplet
	l++;
	
	Asparse->i[l] = pos;
	Asparse->p[l] = pos;
	Asparse->x[l] = 1/temp->value;
      } 
    }
    //take a power source element
    else if(tmp2->type_id == 'i'){
      
      Power_Source* temp = (Power_Source*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	B_array_cs[pos] -= temp->value;	//build the A_array
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	B_array_cs[neg] += temp->value;	//build the A_array
      }
      
    }

    tmp2=tmp2->next;

  }

  C_comp = cs_compress(Asparse);	//create compressed column matrix
    
  cs_dupl(C_comp);
  
      

  //print_A_B_array_cs();
  
  
}


void mna_sparse_tran(int n, int m2, int n_z){
  int pos, neg, k=0, l=-1;
  
  
  size = n + m2;
  
  Csparse = cs_spalloc(size, size, n_z , 1, 1);
  Csparse->nz = n_z;

  element *tmp2 = NULL;

  tmp2 = head;

  while(tmp2!=NULL){

    //take a volt source element
    if(tmp2->type_id== 'v'){
      k++;
    }
    //take an iductance element
    else if(tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(tmp2->type);
      
      //build triplet
      l++;
      
      Csparse->i[l] = n + k;
      Csparse->p[l] = n + k;
      Csparse->x[l] = -temp->value;
      
      k++;
    }
    //take a capacity element
    else if(tmp2->type_id == 'c'){
      
      Capacity* temp = (Capacity*)(tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	//build triplet
	l++;
	
	Csparse->i[l] = pos;
	Csparse->p[l] = pos;
	Csparse->x[l] = temp->value;
	
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	//build triplet
	l++;
	
	Csparse->i[l] = neg;
	Csparse->p[l] = neg;
	Csparse->x[l] = temp->value;
	
	
      }
      
      if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	
	//build triplet
	l++;
	
	Csparse->i[l] = neg;
	Csparse->p[l] = pos;
	Csparse->x[l] = -temp->value;
	
	//build triplet
	l++;
	
	Csparse->i[l] = pos;
	Csparse->p[l] = neg;
	Csparse->x[l] = -temp->value;
	
      }
    
    }
   

    tmp2=tmp2->next;

  }

  C_tilda_comp = cs_compress(Csparse);	//create compressed column matrix
    
  cs_dupl(C_tilda_comp);
 
  //cs_print(Csparse, "C-Sparse.txt", 0);
  
  //cs_print(C_tilda_comp, "C_tilda-compressed-nondupl.txt", 0);
  
  //printf("\nMNA matrix in C_tilda-compressed-nondupl file...\n\n");
  
  
}


void factorization_cs(int flag){
  
  
  
  if(flag == 0){	//lu factorization
    printf("\n\n***** LU method for sparse matrices *****\n\n");
    
    S = cs_sqr(2, C_comp, 0);
    N = cs_lu(C_comp, S, 1);
    cs_spfree(C_comp);
  }
  else if(flag == 1){	//cholesky factorization
    printf("\n\n***** Cholesky method for sparse matrices *****\n\n");
    
    S = cs_schol(1, C_comp);
    N = cs_chol(C_comp, S);
    cs_spfree(C_comp);
  }
}

void lu_solve_cs(int flag){
  //flag gia ton elegxo an exoume dc sweep i oxi, giati i cs_pvec allazei to b dianisma
  if(flag == 0){
    B_array_cs_temp = (double*)calloc(size, sizeof(double));
    
    //initialize x vector
    x_vector = (double*)calloc(size, sizeof(double));
  } 
   
  cs_ipvec(N->pinv, B_array_cs, x_vector, size);
  cs_lsolve(N->L, x_vector);
  cs_usolve(N->U, x_vector);
  cs_ipvec(S->q, x_vector, B_array_cs_temp, size);
    
}

void chol_solve_cs(int flag){
  //flag gia ton elegxo an exoume dc sweep i oxi, giati i cs_pvec allazei to b dianisma
  if(flag == 0){
    B_array_cs_temp = (double*)calloc(size, sizeof(double));
  }
  
  cs_ipvec(S->pinv, B_array_cs, x_vector, size);
  cs_lsolve(N->L, x_vector);
  cs_ltsolve(N->L, x_vector);
  cs_pvec(S->pinv, x_vector, B_array_cs_temp, size);

}

void print_A_B_array_cs(){
  int i;
  
  cs_print(Asparse, "A-Sparse.txt", 0);
  
  cs_print(C_comp, "C-compressed-nondupl.txt", 0);
  
  printf("\nMNA matrix in C-compressed-nondupl file...\n\n");
  
  printf("\nThe %dx1 B array is:\n\n", size);
  
  for(i=0; i < size; i++){
    printf("%.4lf\n", B_array_cs[i]);
  }
  printf("\n\n");
}

void printX_cs(){
  int i;
  //print the x vector
  printf ("x = \n\n");
  
  
  for(i=0; i<size; i++){
    printf("%e\n", B_array_cs_temp[i]);
  }
  
  printf("\n\n");
}