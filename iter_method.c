#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "iter_method.h"
#define MAX_ITERATION 20000
#define EPS 0.00000000000001

double *r,*z,*q,*p, *iter_x, *temp;
double *r_tilda = NULL;
double *z_tilda = NULL;
double *q_tilda = NULL;
double *p_tilda = NULL; 
int size;
int iter;

void cg_method(int cs_flag, int flag, int reset, double itol, int n, int m2){
  double b_norm, r_norm, rho_1, rho, beta, alpha, res;
  int i=0;
  
  iter = 0;
  
  //take the size of the arrays
  size = n + m2;
  
  if(reset == 0){
    //initialize all vectors
    init_cg(cs_flag, flag);
  }
  else{	//keep the last value of x vector for faster convergence
    
    temp = (double*)calloc(size, sizeof(double));
    memcpy(temp, iter_x, size*sizeof(double));
    
    freeIter();
    //initialize all vectors
    init_cg(cs_flag, flag);
    
    memcpy(iter_x, temp, size*sizeof(double));
    
    
    free(temp);
  }
  
  if(cs_flag == 0){
    //Bi-Cg method
    if(flag == 0){
      //r = b - A*iter_x
      dgemv(0, (-1.0), A_array, iter_x, r); 
      vector_add(1.0, r, B_array, r);
    
      //r_tilda = r
      memcpy(r_tilda, r, size*sizeof(double));
      
      //norm of b vector-- only once computed
      b_norm = norm_2(B_array);
      
     
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r vector
      r_norm = norm_2(r);
 
      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
    
	preconditioner_solve(M, r, z);
	
	preconditioner_solve(M, r_tilda, z_tilda);
	
	
	//rho = r^T*z
	rho = ddot(r_tilda, z);
	
	
	
	if(fabs(rho) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
	
	//printf("oook\n");
	if(iter == 1){
	  // p = z
	  memcpy(p, z, (size * sizeof(double)));
	  
	  // p_tilda = z_tilda
	  memcpy(p_tilda, z_tilda, (size * sizeof(double)));
	  
	}
	else{
	  beta = rho / rho_1;
	  
	  //p=z+beta*p
	  vector_add(beta, p, z, p);
	  
	  //p_tilda=z_tilda+beta*p_tilda
	  vector_add(beta, p_tilda, z_tilda, p_tilda);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv(0, 1.0, A_array, p, q); 
	
	//q_tilda = A^T*p_tilda
	dgemv(1, 1.0, A_array, p_tilda, q_tilda); 
	
	//p^T * q  
	res = ddot(p_tilda, q);
	
	//for(i=0; i < size; i++){
	 // printf("%g ... %g \n", q[i], p_tilda[i]);
	//}
	//printf("%g\n", res);
	if(fabs(res) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
      
	alpha = rho / res;
	
	//iter_x = iter_x + alpha*p
	
	vector_add(alpha, p, iter_x, iter_x);
	
	
	//r = r - alpha*q
	vector_add((-alpha), q, r, r);
	
	//r_tilda = r_tilda - alpha*q_tilda
	vector_add((-alpha), q_tilda, r_tilda, r_tilda);
	
	//re-compute the 2nd norm of r vector
	r_norm = norm_2(r);
	
      }
    }
    //Cg method
    else if(flag == 1){
    //r = b - A*iter_x
    dgemv(0, (-1.0), A_array, iter_x, r); 
    vector_add(1.0, r, B_array, r);
  
    
    //norm of b vector-- only once computed
    b_norm = norm_2(B_array);
    
    
    //non zero division condition
    if(b_norm == 0){
      b_norm = 1;
    }
    
    //norm of r vector
    r_norm = norm_2(r);

    while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
      iter++;
      
      preconditioner_solve(M, r, z);
      
      
      //rho = r^T*z
      rho = ddot(r, z);
      
      if(iter == 1){
	// p = z
	memcpy(p, z, (size * sizeof(double)));
	
    
      }
      else{
	beta = rho / rho_1;
	
	//p=z+beta*p
	vector_add(beta, p, z, p);
	
      }
      
      rho_1 = rho;
      
      //q = A*p
      dgemv(0, 1.0, A_array, p, q); 
      
      //p^T * q  
      res = ddot(p, q);
      
    
      alpha = rho / res;
      
      //iter_x = iter_x + alpha*p
      
      vector_add(alpha, p, iter_x, iter_x);
      
      
      //r = r - alpha*q
      vector_add((-alpha), q, r, r);
      
      //re-compute the 2nd norm of r vector
      r_norm = norm_2(r);
      
    }
  }
  }
  else if(cs_flag == 1){
    
      //Bi-Cg method
    if(flag == 0){
      //r = b - A*iter_x
      dgemv_cs(0, (-1.0), C_comp, iter_x, r); 
      vector_add(1.0, r, B_array_cs, r);
    
      //r_tilda = r
      memcpy(r_tilda, r, size*sizeof(double));
      
      //norm of b vector-- only once computed
      b_norm = norm_2(B_array_cs);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r vector
      r_norm = norm_2(r);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
	
	
	preconditioner_solve(M, r, z);
	
	preconditioner_solve(M, r_tilda, z_tilda);
	
	
	//rho = r^T*z
	rho = ddot(r_tilda, z);
	
	if(fabs(rho) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
	
	if(iter == 1){
	  // p = z
	  memcpy(p, z, (size * sizeof(double)));
	  
	  // p_tilda = z_tilda
	  memcpy(p_tilda, z_tilda, (size * sizeof(double)));
	  
	}
	else{
	  beta = rho / rho_1;
	  
	  //p=z+beta*p
	  vector_add(beta, p, z, p);
	  
	  //p_tilda=z_tilda+beta*p_tilda
	  vector_add(beta, p_tilda, z_tilda, p_tilda);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_cs(0, 1.0, C_comp, p, q); 
	
	//q_tilda = A^T*p_tilda
	dgemv_cs(1, 1.0, C_comp, p_tilda, q_tilda); 
	
	//p^T * q  
	res = ddot(p_tilda, q);
	
	if(fabs(res) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
      
	alpha = rho / res;
	
	//iter_x = iter_x + alpha*p
	
	vector_add(alpha, p, iter_x, iter_x);
	
	
	//r = r - alpha*q
	vector_add((-alpha), q, r, r);
	
	//r_tilda = r_tilda - alpha*q_tilda
	vector_add((-alpha), q_tilda, r_tilda, r_tilda);
	
	//re-compute the 2nd norm of r vector
	r_norm = norm_2(r);
	
      }
    }
    //Cg method
    else if(flag == 1){
      //r = b - A*iter_x
      dgemv_cs(0, (-1.0), C_comp, iter_x, r); 
      vector_add(1.0, r, B_array_cs, r);
    
      
      //norm of b vector-- only once computed
      b_norm = norm_2(B_array_cs);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r vector
      r_norm = norm_2(r);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
	
	
	preconditioner_solve(M, r, z);
	
	
	//rho = r^T*z
	rho = ddot(r, z);
	
	if(iter == 1){
	  // p = z
	  memcpy(p, z, (size * sizeof(double)));
	  
      
	}
	else{
	  beta = rho / rho_1;
	  
	  //p=z+beta*p
	  vector_add(beta, p, z, p);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_cs(0, 1.0, C_comp, p, q); 
	
	//p^T * q  
	res = ddot(p, q);
	
      
	alpha = rho / res;
	
	//iter_x = iter_x + alpha*p
	
	vector_add(alpha, p, iter_x, iter_x);
	
	
	//r = r - alpha*q
	vector_add((-alpha), q, r, r);
	
	//re-compute the 2nd norm of r vector
	r_norm = norm_2(r);
	
      }
    }
  }
}

void dgemv(int flag, double a, double *array, double *xx, double *yy){
    
    int i,j;
    if(flag == 0){	//no traspose matrix
      for(i=0; i<size; i++){
	yy[i] = 0;
	for(j=0; j<size; j++){
	  yy[i] += (a * array[(i*size) + j]) * xx[j];
	}
      }
    }
    else if(flag == 1){	//trapsose matrix
      for(i=0; i<size; i++){
	yy[i] = 0;
	for(j=0; j<size; j++){
	  yy[i] += (a * array[(j*size) + i]) * xx[j];
	}
      }
    }
}

void dgemv_cs(int flag, double a, cs *array, double *xx, double *yy){
  int j, pp;
  
  for(j = 0; j<size; j++){
    yy[j] = 0;
  }
  
  if(flag == 0){	//no trapsose matrix
    for (j = 0 ; j < size ; j++){
      for (pp = array->p[j] ; pp < array->p[j+1] ; pp++){
	yy[array->i[pp]] += (a * array->x[pp]) * xx[j] ;
      }
    }
  }
  else if(flag == 1){
    for(j = 0; j<size; j++){
      for(pp = array->p[j] ; pp < array->p[j+1] ; pp++){
	yy[j] += (a * array->x[pp]) * xx[array->i[pp]] ;
      }
    }
  }
}


void vector_add(double a, double *xx, double *yy, double *zz){
  int i;
  
  for(i=0; i<size; i++){
    zz[i] = (a *xx[i]) + yy[i];
  }
}

double norm_2(double *xx){
  double res=0;
  int i;
  
  for(i=0; i<size; i++){
    res += pow(xx[i], 2);
  }
  return sqrt(res);
}

double ddot(double *xx, double *yy){
  double res = 0;
  int i;
  
  for(i=0; i<size; i++){
    res += (xx[i] * yy[i]);
  }
  
  return res;
}

void preconditioner_solve(double *xx, double *yy, double *zz){
  int i;
  
  for(i=0; i<size; i++){
    
    zz[i] = yy[i] / xx[i];
  }
  
}


void init_cg(int cs_flag, int flag){
  int i,j;
  
  //allocate memory for preconditioner M vector
  M = (double*)calloc(size, sizeof(double));
  
  if(M == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  if(cs_flag == 0){
  //M has the diagonal elements of MNA array
    for(i=0; i<size; i++){
      if(A_array[(i*size) + i] == 0){
	M[i] = 1;
      }
      else{
	M[i] =  A_array[(i*size) + i];
      }
    }
  }
  else if(cs_flag == 1){
    int pp;
    
    for(j = 0; j < size; ++j){
      for(pp = C_comp->p[j]; pp < C_comp->p[j+1]; pp++){
	if(C_comp->i[pp] == j){
	  M[j] = C_comp->x[pp];
	  break;
	}
      }	
      if(C_comp->i[pp] != j){
	M[j] = 1; 
      }
    }	
  }
  
  
  //allocate memory for preconditioner r vector
  r = (double*)calloc(size, sizeof(double));
  
  //allocate memory for preconditioner z vector
  z = (double*)calloc(size, sizeof(double));
  
  
  //allocate memory for preconditioner p vector
  p = (double*)calloc(size, sizeof(double));
  
  
  //allocate memory for preconditioner q vector
  q = (double*)calloc(size, sizeof(double));
  
  //allocate memory for preconditioner iter_x vector
  iter_x = (double*)calloc(size, sizeof(double));
  
  if((iter_x == NULL) || (r == NULL) || (z == NULL) || (p == NULL) || (q == NULL)){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  if (flag == 0){
    //allocate memory for preconditioner r_tilda vector
    r_tilda = (double*)calloc(size, sizeof(double));
    
    //allocate memory for preconditioner z_tilda vector
    z_tilda = (double*)calloc(size, sizeof(double));
    
    
    //allocate memory for preconditioner p_tilda vector
    p_tilda = (double*)calloc(size, sizeof(double));
    
    
    //allocate memory for preconditioner q_tilda vector
    q_tilda = (double*)calloc(size, sizeof(double));
    
    if((r_tilda == NULL) || (z_tilda == NULL) || (p_tilda == NULL) || (q_tilda == NULL)){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
    }
  }
  
}

void printIter_x(int flag){
  int i;
  
  if(flag == 0){
    
    printf("\n*****Bi-CG METHOD*****\n\n");
    printf("The Bi-CG algorith ran %d times until solution...\n\n", iter);
    //print the x vector
    printf ("x = \n\n");
    
    for(i=0; i<size; i++){
      printf("%e\n", iter_x[i]);
    }
    
    printf("\n\n");
  }
  else if(flag == 1){
    printf("\n*****CG METHOD*****\n\n");
    printf("The CG algorith ran %d times until solution...\n\n", iter);
    //print the x vector
    printf ("x = \n\n");
    
    for(i=0; i<size; i++){
      printf("%e\n", iter_x[i]);
    }
    
    printf("\n\n");
  }
}


void freeIter(){
  
  free(iter_x);
  
  free(r);
  
  free(z);
  
  free(q);
  
  free(p);
  
  free(M);
  
  free(r_tilda);
  
  free(z_tilda);
  
  free(q_tilda);
  
  free(p_tilda);
  
    
}