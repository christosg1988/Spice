#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "dc_sweep.h"


int sweep_size;
int pos = -1;
int neg = -1;


gsl_vector *temp_x = NULL;
gsl_permutation * temp_p = NULL;
double *x_cs = NULL;
double *temp2_array = NULL;

void dc_sweep(double start_value, double end_value, double inc_value){
  double i;
  int j, m, l;
  //double g;
  j = 0;
  
  temp_x = x;
  temp_p = p;
  temp2_array = tmp2_array;
  
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
  
  
  init_dc_sweep_array();
  
  //LU decomposition for dc sweep
  
  i=start_value;  
  for(l=0; l < sweep_size; l++){    
    temp2_array[n + sweep_k] = i;
   
    lu_solve();
    
    //keep the first n elements of x vector
    for(m=0; m < size; m++){
      
      dc_sweep_array[j] = gsl_vector_get(temp_x, m);
      j++;
    }
    printf("\n");
    i += inc_value;
  }
}


void dc_sweep_cs(double start_value, double end_value, double inc_value){
  double i;
  int j, m, l;
  //double g;
  j = 0;
  
  //x_cs = x_vector;
  //temp2_array = B_array_cs;
  
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
 
  init_dc_sweep_array();
  
  //LU decomposition for dc sweep
  
  i=start_value;  
  for(l=0; l < sweep_size; l++){ 
  
    B_array_cs[n + sweep_k] = i;
    
    
    lu_solve_cs(1);
    
    //keep the first n elements of x vector
    for(m=0; m < size; m++){
      dc_sweep_array[j] = B_array_cs_temp[m];
      j++;
    }
    
  
    i += inc_value;
  }
}


void iter_dc_sweep(char t, int flag, int n, int m2, double start_value, double end_value, double inc_value){
  double i, itol;
  int m, l;
  int j = 0;
  
  itol = itol_t;
   
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
  
  init_dc_sweep_array();
  
  
  i=start_value;  
  for(l=0; l < sweep_size; l++){  
    
    if(t == 'v'){
      B_array[n + sweep_k] = i;
    }
    else if(t == 'i'){
      if(pos != -1){
	B_array[pos] = -i;
      }
    
      //check if negative node is not the ground
      if(neg != -1){
	B_array[neg] = i;
      }
    }
    
    cg_method(0, flag, j, itol, n, m2);
    
    //keep the first n elements of x vector
    for(m=0; m < size; m++){
      dc_sweep_array[j] = iter_x[m];
      j++;
    }
    
    i += inc_value;
    
  }
  
}


void iter_dc_sweep_cs(char t, int flag, int n, int m2, double start_value, double end_value, double inc_value){
  double i, itol;
  int m, l;
  int j = 0;
  
  itol = itol_t;
   
  
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
  
  init_dc_sweep_array();
  
  
  i=start_value;  
  for(l=0; l < sweep_size; l++){  
    
    if(t == 'v'){
      B_array_cs[n + sweep_k] = i;
    }
    else if(t == 'i'){
      if(pos != -1){
	B_array_cs[pos] = -i;
      }
    
      //check if negative node is not the ground
      if(neg != -1){
	B_array_cs[neg] = i;
      }
    }
    
    cg_method(1, flag, j, itol, n, m2);
    
    
    //keep the first n elements of x vector
    for(m=0; m < size; m++){
      dc_sweep_array[j] = iter_x[m];
      j++;
    }
    
    i += inc_value;
    
  }
  
}


void dc_sweep2(double start_value, double end_value, double inc_value, int flag){
  double i;
  int j, m, l;
  //double g;
  j = 0;
  
  temp_x = x;
  temp_p = p;
  temp2_array = tmp2_array;
  
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
  
  init_dc_sweep_array();
  
  i=start_value; 
  if(flag == 0){// LU decomposition for dc sweep
    for(l=0; l < sweep_size; l++){
      //check if positive node is not the ground
      if(pos != -1){
	temp2_array[pos] = -i;
      }
      
      //check if negative node is not the ground
      if(neg != -1){
	temp2_array[neg] = i;
      }
      
      lu_solve();
      
      //keep the first n elements of x vector
      for(m=0; m < size; m++){
	dc_sweep_array[j] = gsl_vector_get(temp_x, m);
	j++;
      }
      i += inc_value;
    }
  }
  else if(flag == 1){//Cholesky method for dc sweep
    for(l=0; l < sweep_size; l++){
      
      //check if positive node is not the ground
      if(pos != -1){
	temp2_array[pos] = -i;
      }
      
      //check if negative node is not the ground
      if(neg != -1){
	temp2_array[neg] = i;
      }
      
      chol_solve();
      
      //keep the first n elements of x vector
      for(m=0; m < size; m++){
	dc_sweep_array[j] = gsl_vector_get(temp_x, m);
	j++;
      }
      i += inc_value;
    }
  }
}


void dc_sweep2_cs(double start_value, double end_value, double inc_value, int flag){
  double i;
  int j, m, l;
  //double g;
  j = 0;

  
  //take the number of for-loops for the dc weep
  sweep_size = (int)(((end_value - start_value) / inc_value) + 1.0);
  
  init_dc_sweep_array();
  
  i=start_value; 
  if(flag == 0){// LU decomposition for dc sweep
    for(l=0; l < sweep_size; l++){
      //check if positive node is not the ground
      if(pos != -1){
	B_array_cs[pos] = -i;
      }
      
      //check if negative node is not the ground
      if(neg != -1){
	B_array_cs[neg] = i;
      }
      
      lu_solve_cs(1);
      
      //keep the first n elements of x vector
      for(m=0; m < size; m++){
	dc_sweep_array[j] =  B_array_cs_temp[m];
	j++;
      }
      i += inc_value;
    }
  }
  else if(flag == 1){//Cholesky method for dc sweep
    for(l=0; l < sweep_size; l++){
      
      //check if positive node is not the ground
      if(pos != -1){
	B_array_cs[pos] = -i;
      }
      
      //check if negative node is not the ground
      if(neg != -1){
	B_array_cs[neg] = i;
      }
      
      chol_solve_cs(1);
      
      //keep the first n elements of x vector
      for(m=0; m < size; m++){
	dc_sweep_array[j] = B_array_cs_temp[m];
	j++;
      }
      i += inc_value;
    }
  }
}


void dc_sweep_printf(char *str, char c){
  int i, j;
  
  if(c == 'v'){
    printf("DC sweep for V%s:\n\n", str);
  }
  else if(c == 'i'){
    printf("DC sweep for I%s:\n\n", str);
  }
  
  //print all the x vectors for the dc sweep
  for(i=0; i < sweep_size; i++){
    for(j=0; j < size; j++){
      printf("x%d_%d = %.5lf\t", (i+1), (j+1), dc_sweep_array[(i*size) + j]);
    }
    printf("\n");
  }
  
  printf("\n\n");
}

void init_dc_sweep_array(){
  //initialize the output array for dc sweep
  dc_sweep_array = (double*)calloc(sweep_size*size, sizeof(double));
}

void free_dc_sweep_array(){
  free(dc_sweep_array);
}

void search_list(char *str){
  int k=0;
  element *src = head;
  
  while(src != NULL){
    
    if(src->type_id== 'v'){
      Volt_Source* temp = (Volt_Source*)(src->type);
  
      if(strcmp(str, temp->name) == 0){
	sweep_k = k;
	break;
      }
      k++;
    }
    
    if(src->type_id== 'l'){	//gia na vrei tin sosti thesi tis tasis sto mna auksanei to k (stoixeia tis omadas 2), akoma kai an sinantisei L
      k++;
    }
    
    src = src->next;
  }
  
}

void search_list_nodes(char *str){
  element *src2 = head;
  
  while(src2 != NULL){
    
    if(src2->type_id== 'i'){
      Power_Source* temp = (Power_Source*)(src2->type);
  
      if(strcmp(str, temp->name) == 0){
	if(strcmp(temp->pos_node, "0") != 0){
	  pos = search_element_id(temp->pos_node);
	}
	if(strcmp(temp->neg_node, "0") != 0){
	  neg = search_element_id(temp->neg_node);
	}
	break;
      }
    }
    
    src2 = src2->next;
  }
}