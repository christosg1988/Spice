#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "parser.h"
#include "hashtable.h"
#include "mna.h"
#include "factorization.h"


gsl_vector *x = NULL;
gsl_permutation * p = NULL;
double *tmp_array = NULL;
double *tmp2_array = NULL;

int size = 0;

void factorization(int flag, int n, int m2){
  //take the total size of the array
  size = n+m2;
 
  //flag for SPD option in netlist
  if(flag == 0){
    
    init_LU();
    
    lu_Decomp();
    
    printf("*****LU decomposition*****\n\n");
  }
  else if(flag == 1){
    
    init_Chol();
    chol_Decomp();
    
    printf("*****Cholesky decomposition*****\n\n");
  }
}

void init_LU(){
  int i;
  
  //allocate memory for x vector
  x = gsl_vector_calloc(size);
  
  //allocate memory  for p permutatation
  p = gsl_permutation_calloc (size);
  
  //allocate memory for a temp A array 
  tmp_array = (double*)calloc((size*size), sizeof(double));
  
  //copy all elements of A array to temp array
  for(i=0; i < size*size; i++){
    tmp_array[i] = A_array[i];
  }
  
  //allocate memory for a temp b vector
  tmp2_array = (double*)calloc(size, sizeof(double));
  
  //copy all elements of b vector to temp2 array
  for(i=0; i < size; i++){
    tmp2_array[i] = B_array[i];
  }
}

void init_Chol(){
  int i;
  
  //allocate memory for x vector
  x = gsl_vector_calloc(size);
  
  //allocate memory for a temp A array 
  tmp_array = (double*)calloc((size*size), sizeof(double));
  
  //copy all elements of A array to temp array
  for(i=0; i < size*size; i++){
    tmp_array[i] = A_array[i];
  }
  
  //allocate memory for a temp b vector
  tmp2_array = (double*)calloc(size, sizeof(double));
  
  //copy all elements of b vector to temp2 array
  for(i=0; i < size; i++){
    tmp2_array[i] = B_array[i];
  }
}

void lu_Decomp(){
  
  int s;
  
  //take the temp vector as gsl view matrix into scope
  gsl_matrix_view m = gsl_matrix_view_array(tmp_array, size, size);
  
  
  //LU decomposition of mna (also permutatation vector calculated)
  gsl_linalg_LU_decomp (&m.matrix, p, &s);

 
}

void lu_solve(){
  //take the b vector as gsl view vector into scope
  gsl_vector_view b = gsl_vector_view_array(tmp2_array, size);
  
  //take the temp vector as gsl view matrix into scope
  gsl_matrix_view m = gsl_matrix_view_array(tmp_array, size, size);
  
  gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
}

void chol_Decomp(){
  
  //take the temp vector as gsl view matrix into scope
  gsl_matrix_view m = gsl_matrix_view_array(tmp_array, size, size);
  
  gsl_linalg_cholesky_decomp(&m.matrix);
  
}

void chol_solve(){
  //take the b vector as gsl view vector into scope
  gsl_vector_view b = gsl_vector_view_array(tmp2_array, size);
  
  //take the temp vector as gsl view matrix into scope
  gsl_matrix_view m = gsl_matrix_view_array(tmp_array, size, size);
  
  gsl_linalg_cholesky_solve(&m.matrix, &b.vector, x);
  
}

void printChol(){
  int i, j;
  
  
  
  //print the L array
  printf("L = \n\n");
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      if(i > j){
	printf("0.0000\t");
      }
      else{
	printf("%.4lf\t", tmp_array[(j*size) + i]);
      }
    }
    printf("\n");
  }
  
  printf("\n\n");
  
  //print the L^T array
  printf("L^T = \n\n");
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      if(j > i){
	printf("0.0000\t");
      }
      else{
	printf("%.4lf\t", tmp_array[(j*size) + i]);
      }
    }
    printf("\n");
  }
  
  printf("\n\n");
  
}


void printLU(){
  int i, j;
  
 
  
  //print the L array
  printf("L = \n\n");
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      if(j == i){
	printf("1.0000\t");
      }
      else if(i > j){
	printf("0.0000\t");
      }
      else{
	printf("%.4lf\t", tmp_array[(j*size) + i]);
      }
    }
    printf("\n");
  }
  
  printf("\n\n");
  
  //print the U array
  printf("U = \n\n");
  
  for(j=0; j < size; j++){
    for(i=0; i < size; i++){
      if(j > i){
	printf("0.0000\t");
      }
      else{
	printf("%.4lf\t", tmp_array[(j*size) + i]);
      }
    }
    printf("\n");
  }
  
  printf("\n\n");
  
  
}

void printX(){
  
  //print the x vector
  printf ("x = \n\n");
  
  gsl_vector_fprintf(stdout, x, "%e");
  
  printf("\n\n");
}

void freeLU(){
  
  gsl_vector_free (x);
  
  gsl_permutation_free (p);
  
  free(tmp_array);
  
  free(tmp2_array);
}

