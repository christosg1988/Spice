#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex_math.h>
#include "parser.h"
#include "hashtable.h"
#include "ac_analysis.h"
#define PI 3.14159265

element *ac_tmp2 = NULL;
int size;

gsl_matrix_complex *A_ac = NULL;
gsl_vector_complex *x_ac = NULL;
gsl_permutation *p_ac = NULL;
gsl_vector_complex *b_ac = NULL;

void init_A(){
  //initialize two-dimensional A array
  A_ac_array = (complex double*)calloc((size*size), sizeof(complex double));
}

void init_A_sparse(int n_z){
  //initialize two-dimensional A array
  A_ac_sparse = cs_ci_spalloc(size, size, n_z , 1, 1);
  A_ac_sparse->nz = n_z;
}


void init_B(){
  //initialize one-dimensional B array
  B_ac_array = (complex double*)calloc(size, sizeof(complex double));
}

/*
void mna_ac_cs(int n, int m2, double freq, int n_z){
  int pos, neg, k=0, l=-1;
  
  size = n +m2;
  
  //initialize
  Asparse_ac = cs_cl_spalloc(size, size, n_z, 1, 1);
}

*/
void mna_ac(int n, int m2, double freq){
  int k = 0;
  int pos, neg;
  complex double rad_freq=0;
  
  size = n + m2;	//take the size of the A and B array
  
  rad_freq = 2 * PI * freq;
  rad_freq *= I;	//take the result j*ω
  
  ac_tmp2 = head;
  
  while(ac_tmp2!=NULL){
    
    //take a volt source element
    if(ac_tmp2->type_id== 'v'){
      
      Volt_Source* temp = (Volt_Source*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	A_ac_array[((n + k)*size) + pos] += 1;	//build the A_ac_array
	A_ac_array[(pos*size) + (n + k)] += 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	A_ac_array[((n + k)*size) + neg] -= 1;	//build the A_ac_array
	A_ac_array[(neg*size) + (n + k)] -= 1;
      }
      
      B_ac_array[n+k] += (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the B_ac_array
      
      k++;
    
    }
    //take an iductance element
    else if(ac_tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	A_ac_array[(n + k)*size + pos] += 1;	//build the A_ac_array
	A_ac_array[(pos*size) + (n + k)] += 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	A_ac_array[(n + k)*size + neg] -= 1;	//build the A_ac_array
	A_ac_array[(neg*size) + (n + k)] -= 1;
      }
      
      A_ac_array[(n + k)*size + (n + k)] -= rad_freq * temp->value;
      k++;
    }
    //take a resistance element
    else if(ac_tmp2->type_id == 'r'){
      
      Resistance* temp = (Resistance*)(ac_tmp2->type);
      
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
	A_ac_array[(pos*size) + neg] -= (1/(temp->value));		//build the A_ac_array
	A_ac_array[(neg*size) + pos] -= (1/(temp->value));
	A_ac_array[(pos*size) + pos] += (1/(temp->value));
	A_ac_array[(neg*size) + neg] += (1/(temp->value));
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	A_ac_array[(neg*size) + neg] += (1/(temp->value));
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	A_ac_array[(pos*size) + pos] += (1/(temp->value));
      } 
    }
    //take a power source element
    else if(ac_tmp2->type_id == 'i'){
      
      Power_Source* temp = (Power_Source*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	//printf("re\n");
	B_ac_array[pos] -= (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the A_ac_array
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	B_ac_array[neg] += (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the A_ac_array
	//printf("%g\n",temp->value);
      }
      
    }
    else if(ac_tmp2->type_id == 'c'){
      Capacity* temp = (Capacity*)(ac_tmp2->type);
      
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
	A_ac_array[(pos*size) + neg] -= rad_freq * temp->value;		//build the A_ac_array
	A_ac_array[(neg*size) + pos] -= rad_freq * temp->value;
	A_ac_array[(pos*size) + pos] += rad_freq * temp->value;
	A_ac_array[(neg*size) + neg] += rad_freq * temp->value;
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	A_ac_array[(neg*size) + neg] += rad_freq * temp->value;
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	A_ac_array[(pos*size) + pos] += rad_freq * temp->value;
      } 
    }
    ac_tmp2=ac_tmp2->next;
  
  }
  
}



void mna_ac_sparse(int n, int m2, double freq){
  int l = -1, k=0;
  int pos, neg;
  complex double rad_freq=0;
  
  size = n + m2;	//take the size of the A and B array
  
  rad_freq = 2 * PI * freq;
  rad_freq *= I;	//take the result j*ω
  
  ac_tmp2 = head;
  
  while(ac_tmp2!=NULL){
    
    //take a volt source element
    if(ac_tmp2->type_id== 'v'){
      
      Volt_Source* temp = (Volt_Source*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = n + k;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = 1;
	
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = n + k;
	A_ac_sparse->x[l] = 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = n + k;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = -1;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = n + k;
	A_ac_sparse->x[l] = -1;
      }
      
      B_ac_array[n+k] += (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the B_ac_array
      
      k++;
    
    }
    //take an iductance element
    else if(ac_tmp2->type_id== 'l'){
      
      Inductance* temp = (Inductance*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = n + k;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = 1;
	
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = n + k;
	A_ac_sparse->x[l] = 1;
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = n + k;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = -1;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = n + k;
	A_ac_sparse->x[l] = -1;
      }
      l++;
      
      A_ac_sparse->i[l] = n + k;
      A_ac_sparse->p[l] = n + k;
      A_ac_sparse->x[l] = -(rad_freq * temp->value);
      
      k++;
    }
    //take a resistance element
    else if(ac_tmp2->type_id == 'r'){
      
      Resistance* temp = (Resistance*)(ac_tmp2->type);
      
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
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = -1/temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = -1/temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = 1/temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = 1/temp->value;
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = 1/temp->value;
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	//build triplet
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = 1/temp->value;
      } 
    }
    //take a power source element
    else if(ac_tmp2->type_id == 'i'){
      
      Power_Source* temp = (Power_Source*)(ac_tmp2->type);
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	//printf("re\n");
	B_ac_array[pos] -= (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the A_ac_array
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	B_ac_array[neg] += (temp->ac_ptr->mag * (cos(2 * PI * temp->ac_ptr->phase/360) + I * sin(2 * PI * temp->ac_ptr->phase/360)));	//build the A_ac_array
	//printf("%g\n",temp->value);
      }
      
    }
    else if(ac_tmp2->type_id == 'c'){
      Capacity* temp = (Capacity*)(ac_tmp2->type);
      
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
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = -rad_freq * temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = -rad_freq * temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = rad_freq * temp->value;
	
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = rad_freq * temp->value;
	
      }
      //in case that positive node is the ground node
      else if(strcmp(temp->pos_node, "0") == 0){
	l++;
	
	A_ac_sparse->i[l] = neg;
	A_ac_sparse->p[l] = neg;
	A_ac_sparse->x[l] = rad_freq * temp->value;
      }
      //in case that the negative node is the ground node
      else if(strcmp(temp->neg_node, "0") == 0){
	l++;
	
	A_ac_sparse->i[l] = pos;
	A_ac_sparse->p[l] = pos;
	A_ac_sparse->x[l] = rad_freq * temp->value;
      } 
    }
    ac_tmp2=ac_tmp2->next;
  
  }
  
  C_ac_comp = cs_ci_compress(A_ac_sparse);	//create compressed column matrix
    
  cs_ci_dupl(C_ac_comp);
}


void factorization_cs_ac(int flag){
  
  if(flag == 0){	//lu factorization
    
    S_ac = cs_ci_sqr(2, C_ac_comp, 0);
    N_ac = cs_ci_lu(C_ac_comp, S_ac, 1);
    cs_ci_spfree(C_ac_comp);
  }
  else if(flag == 1){	//cholesky factorization
    
    S_ac = cs_ci_schol(1, C_ac_comp);
    N_ac = cs_ci_chol(C_ac_comp, S_ac);
    cs_ci_spfree(C_ac_comp);
  }
  
  x_ac_cs = (complex double*)calloc(size, sizeof(complex double));
}

void lu_solve_cs_ac(int flag){
  if(flag == 0){
    B_ac_array_temp = (complex double*)calloc(size, sizeof(complex double));
    
    //initialize x vector
    //x_vector = (double*)calloc(size, sizeof(double));
  } 
   
  cs_ci_ipvec(N_ac->pinv, B_ac_array, x_ac_cs, size);
  cs_ci_lsolve(N_ac->L, x_ac_cs);
  cs_ci_usolve(N_ac->U, x_ac_cs);
  cs_ci_ipvec(S_ac->q, x_ac_cs, B_ac_array_temp, size);
    
}

void print_A_B_array_cs_ac(){
  int i;
  
  cs_ci_print(A_ac_sparse, 0);
  
  cs_ci_print(C_ac_comp,  0);
  
  printf("\nMNA matrix in C-compressed-nondupl file...\n\n");
  
  printf("\nThe %dx1 B array is:\n\n", size);
  
  for(i=0; i < size; i++){
    printf("%.4lf, %.4lf\n", creal(B_ac_array[i]), cimag(B_ac_array[i]));
  }
  printf("\n\n");
}

void print_A_ac_array(){
  int i, j;
  
  printf("\nThe %dx%d A array for AC analysis is :\n\n", size, size);
  
  for(j=0; j < size; j++){
    printf("%d:", j);
    for(i=0; i < size; i++){
      printf(" %.4lf + j*(%.4lf)\t", creal(A_ac_array[(j*size) + i]), cimag(A_ac_array[(j*size) + i]));
    }
    printf("\n");
  }
  printf("\n\n");
}

void print_B_ac_array(){
  int i;
  
  printf("\nThe %dx1 B array for AC analysis is:\n\n", size);
  
  for(i=0; i < size; i++){
    printf("%.4lf + j*(%.4lf)\n", creal(B_ac_array[i]), cimag(B_ac_array[i]));
  }
  printf("\n\n");
}

void lu_fact(){
  int i,j,s;
  gsl_complex value;
  
  //in every step the MNA matrix is different
  //copy all elements from double *A to gsl_matrix A
  for(j=0; j<size; j++){
    for(i=0; i<size; i++){
      value = gsl_complex_rect(creal(A_ac_array[(j*size) + i]), cimag(A_ac_array[(j*size) + i]));
      gsl_matrix_complex_set(A_ac, i, j, value);
    }
  }
  
  //LU decomposition of mna (also permutatation vector calculated)
  gsl_linalg_complex_LU_decomp (A_ac, p_ac, &s);

  gsl_linalg_complex_LU_solve (A_ac, p_ac, b_ac, x_ac);
 
}


void chol_fact(){
  int i,j;
  gsl_complex value;
  
  //in every step the MNA matrix is different
  //copy all elements from double *A to gsl_matrix A
  for(j=0; j<size; j++){
    for(i=0; i<size; i++){
      value = gsl_complex_rect(creal(A_ac_array[(j*size) + i]), cimag(A_ac_array[(j*size) + i]));
      gsl_matrix_complex_set(A_ac, i, j, value);
    }
  }
  
  //LU decomposition of mna 
  gsl_linalg_complex_cholesky_decomp (A_ac);

  gsl_linalg_complex_cholesky_solve (A_ac, b_ac, x_ac);
  
 
}



void initialize_lu(){
  int i;
  gsl_complex value;
  //allocate memory for x vector
  x_ac = gsl_vector_complex_calloc(size);
  
  //allocate memory  for p permutatation
  p_ac = gsl_permutation_calloc (size);
  
  A_ac = gsl_matrix_complex_calloc(size, size);
  
  b_ac = gsl_vector_complex_calloc(size);
  
 
  //copy all elements from double *b to gsl_vector b
  for(i=0; i<size; i++){
    value = gsl_complex_rect(creal(B_ac_array[i]), cimag(B_ac_array[i]));
    gsl_vector_complex_set(b_ac, i, value);
  }
 
  
}


void initialize_chol(){
  int i;
  gsl_complex value;
  //allocate memory for x vector
  x_ac = gsl_vector_complex_calloc(size);
  
  A_ac = gsl_matrix_complex_calloc(size, size);
  
  b_ac = gsl_vector_complex_calloc(size);
  
 
  //copy all elements from double *b to gsl_vector b
  for(i=0; i<size; i++){
    value = gsl_complex_rect(creal(B_ac_array[i]), cimag(B_ac_array[i]));
    gsl_vector_complex_set(b_ac, i, value);
  }
 
  
}


void printX_ac(){
  int i;
  //print the x vector
  printf ("x = \n\n");
  for(i=0; i<size; i++){
    printf("%e + j*(%e)\n", GSL_REAL(gsl_vector_complex_get(x_ac, i)), GSL_IMAG(gsl_vector_complex_get(x_ac, i)));
  }
  printf("\n\n");
}

void printX_ac_cs(){
  int i;
  //print the x vector
  printf ("x = \n\n");
  for(i=0; i<size; i++){
    printf("%e + j*(%e)\n", creal(B_ac_array_temp[i]), cimag(B_ac_array_temp[i]));
  }
  printf("\n\n");
}

void freeLU_ac(){
  
  gsl_vector_complex_free (x_ac);
  
  gsl_permutation_free (p_ac);
  
}

void freeChol_ac(){
  
  gsl_vector_complex_free (x_ac);

}
