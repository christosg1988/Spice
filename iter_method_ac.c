#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "iter_method_ac.h"

#define MAX_ITERATION 20000
#define EPS 0.00000000000001

int w = 0;

complex double *r_ac,*z_ac,*q_ac,*p_ac, *iter_x_ac, *temp_ac;
complex double *r_tilda_ac = NULL;
complex double *z_tilda_ac = NULL;
complex double *q_tilda_ac = NULL;
complex double *p_tilda_ac = NULL; 
int size;
int iter;

void cg_method_ac(int cs_flag, int flag, int reset, double itol, int n, int m2){
  double b_norm, r_norm;
  complex double rho_1, rho, beta, alpha, res;
  int i=0;
  
  iter = 0;
  
  //take the size of the arrays
  size = n + m2;
  
  if(reset == 0){
    //initialize all vectors
    init_cg_ac(cs_flag, flag);
  }
  else{	//keep the last value of x vector for faster convergence
    temp_ac = (complex double*)calloc(size, sizeof(complex double));
    memcpy(temp_ac, iter_x_ac, size*sizeof(complex double));
    
    freeIter_ac();
    //initialize all vectors
    init_cg_ac(cs_flag, flag);
    
    memcpy(iter_x_ac, temp_ac, size*sizeof(complex double));
    
    free(temp_ac);
  }
  
  if(cs_flag == 0){
    /*for(i=0; i <size; ++i){
      printf("%g, %g\n", creal(M_ac[i]), cimag(M_ac[i]));
    }*/
    //Bi-Cg method
    if(flag == 0){
      //r_ac = b - A*iter_x_ac
      dgemv_ac(0, (-1.0), A_ac_array, iter_x_ac, r_ac); 
      vector_add_ac(0, 1.0, r_ac, B_ac_array, r_ac);
    
      //r_tilda_ac = r_ac
      memcpy(r_tilda_ac, r_ac, size*sizeof(complex double));
      
      //norm of b vector-- only once computed
      b_norm = norm_2_ac(B_ac_array);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r_ac vector
      r_norm = norm_2_ac(r_ac);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
    
	preconditioner_solve_ac(0, M_ac, r_ac, z_ac);
	
	preconditioner_solve_ac(1, M_ac, r_tilda_ac, z_tilda_ac);
	//printf("ok\n");
	
	//rho = r_ac^H*z_ac
	rho = ddot_ac(r_tilda_ac, z_ac);
	
	//for(i=0; i<size; i++){
	 // printf("%g,%g  ... %g,%g  \n", creal(z_ac[i]), cimag(z_ac[i]), creal(r_tilda_ac[i]), cimag(r_tilda_ac[i]));
       // }
	
	if(fabs(rho) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
	
	if(iter == 1){
	  // p = z_ac
	  memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	  
	  // p_tilda_ac = z_tilda_ac
	  memcpy(p_tilda_ac, z_tilda_ac, (size * sizeof(complex double)));
	  
	}
	else{
	  beta = divide(rho, rho_1);
	  
	  //p=z_ac+beta*p_ac
	  vector_add_ac(0, beta, p_ac, z_ac, p_ac);
	  
	  //p_tilda_ac=z_tilda_ac+beta*p_tilda_ac
	  vector_add_ac(0, conj(beta), p_tilda_ac, z_tilda_ac, p_tilda_ac);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_ac(0, 1.0, A_ac_array, p_ac, q_ac); 
	
	//q_tilda_ac = A^H*p_tilda_ac
	dgemv_ac(1, 1.0, A_ac_array, p_tilda_ac, q_tilda_ac); 
	
	//p^H * q  
	res = ddot_ac(p_tilda_ac, q_ac);
	
	if(fabs(res) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
      
	alpha = divide(rho ,res);
	
	//iter_x_ac = iter_x_ac + alpha*p
	
	vector_add_ac(0, alpha, p_ac, iter_x_ac, iter_x_ac);
	
	
	//r_ac = r_ac - alpha*q
	vector_add_ac(0, -alpha, q_ac, r_ac, r_ac);
	
	//r_tilda_ac = r_tilda_ac - alpha*q_tilda_ac
	vector_add_ac(0, conj(-alpha), q_tilda_ac, r_tilda_ac, r_tilda_ac);
	
	//re-compute the 2nd norm of r_ac vector
	r_norm = norm_2_ac(r_ac);
	
      }
    }
    //Cg method
    else if(flag == 1){
      
    //r_ac = b - A*iter_x_ac
    dgemv_ac(0, (-1.0), A_ac_array, iter_x_ac, r_ac); 
    
    vector_add_ac(0, 1.0, r_ac, B_ac_array, r_ac);
 
  
    //norm of b vector-- only once computed
    b_norm = norm_2_ac(B_ac_array);
     
    
    //non zero division condition
    if(b_norm == 0){
      b_norm = 1;
    }
    
    //norm of r_ac vector
    r_norm = norm_2_ac(r_ac);

    while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
      iter++;
      
      preconditioner_solve_ac(0, M_ac, r_ac, z_ac);
     
      //rho = r_ac^T*z_ac
      rho = ddot_ac(r_ac, z_ac);
       
      
      if(iter == 1){
	// p = z_ac
	memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	
      }
      else{
	beta = divide(rho,rho_1);
	
	//p=z_ac+beta*p
	
	
	vector_add_ac(0, beta, p_ac, z_ac, p_ac);
	
      }
      
      rho_1 = rho;
      
      //q = A*p
      dgemv_ac(0, 1.0, A_ac_array, p_ac, q_ac); 
      
     
     // printf("\n");
      //p^T * q  
      res = ddot_ac(p_ac, q_ac);
      
   
      alpha = divide(rho, res);
      
      //iter_x_ac = iter_x_ac + alpha*p
      vector_add_ac(0, alpha, p_ac, iter_x_ac, iter_x_ac);
     
      //r = r_ac - alpha*q
      vector_add_ac(0, (-alpha), q_ac, r_ac, r_ac);
      
      //re-compute the 2nd norm of r vector
      r_norm = norm_2_ac(r_ac);
      //printf("edw\n");
    }
  }
  }
  
  // Sparse matrices
  if(cs_flag == 1){
    //Bi-Cg method
    if(flag == 0){
      //r_ac = b - A*iter_x_ac
      dgemv_cs_ac(0, (-1.0), C_ac_comp, iter_x_ac, r_ac); 
      vector_add_ac(0, 1.0, r_ac, B_ac_array, r_ac);
    
      //r_tilda_ac = r_ac
      memcpy(r_tilda_ac, r_ac, size*sizeof(complex double));
      
      //norm of b vector-- only once computed
      b_norm = norm_2_ac(B_ac_array);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r_ac vector
      r_norm = norm_2_ac(r_ac);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
    
	preconditioner_solve_ac(0, M_ac, r_ac, z_ac);
	
	preconditioner_solve_ac(1, M_ac, r_tilda_ac, z_tilda_ac);
	//printf("ok\n");
	
	//rho = r_ac^H*z_ac
	rho = ddot_ac(r_tilda_ac, z_ac);
	
	//for(i=0; i<size; i++){
	 // printf("%g,%g  ... %g,%g  \n", creal(z_ac[i]), cimag(z_ac[i]), creal(r_tilda_ac[i]), cimag(r_tilda_ac[i]));
       // }
	
	if(fabs(rho) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
	
	if(iter == 1){
	  // p = z_ac
	  memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	  
	  // p_tilda_ac = z_tilda_ac
	  memcpy(p_tilda_ac, z_tilda_ac, (size * sizeof(complex double)));
	  
	}
	else{
	  beta = divide(rho, rho_1);
	  
	  //p=z_ac+beta*p_ac
	  vector_add_ac(0, beta, p_ac, z_ac, p_ac);
	  
	  //p_tilda_ac=z_tilda_ac+beta*p_tilda_ac
	  vector_add_ac(0, conj(beta), p_tilda_ac, z_tilda_ac, p_tilda_ac);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_cs_ac(0, 1.0, C_ac_comp, p_ac, q_ac); 
	
	//q_tilda_ac = A^H*p_tilda_ac
	dgemv_cs_ac(1, 1.0, C_ac_comp, p_tilda_ac, q_tilda_ac); 
	
	//p^H * q  
	res = ddot_ac(p_tilda_ac, q_ac);
	
	if(fabs(res) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
      
	alpha = divide(rho ,res);
	
	//iter_x_ac = iter_x_ac + alpha*p
	
	vector_add_ac(0, alpha, p_ac, iter_x_ac, iter_x_ac);
	
	
	//r_ac = r_ac - alpha*q
	vector_add_ac(0, -alpha, q_ac, r_ac, r_ac);
	
	//r_tilda_ac = r_tilda_ac - alpha*q_tilda_ac
	vector_add_ac(0, conj(-alpha), q_tilda_ac, r_tilda_ac, r_tilda_ac);
	
	//re-compute the 2nd norm of r_ac vector
	r_norm = norm_2_ac(r_ac);
	
      }
    }
    //Cg method
    else if(flag == 1){
      
    //r_ac = b - A*iter_x_ac
    dgemv_ac(0, (-1.0), A_ac_array, iter_x_ac, r_ac); 
    
    vector_add_ac(0, 1.0, r_ac, B_ac_array, r_ac);
 
  
    //norm of b vector-- only once computed
    b_norm = norm_2_ac(B_ac_array);
     
    
    //non zero division condition
    if(b_norm == 0){
      b_norm = 1;
    }
    
    //norm of r_ac vector
    r_norm = norm_2_ac(r_ac);

    while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
      iter++;
      
      preconditioner_solve_ac(0, M_ac, r_ac, z_ac);
     
      //rho = r_ac^T*z_ac
      rho = ddot_ac(r_ac, z_ac);
       
      
      if(iter == 1){
	// p = z_ac
	memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	
      }
      else{
	beta = divide(rho,rho_1);
	
	//p=z_ac+beta*p
	
	
	vector_add_ac(0, beta, p_ac, z_ac, p_ac);
	
      }
      
      rho_1 = rho;
      
      //q = A*p
      dgemv_ac(0, 1.0, A_ac_array, p_ac, q_ac); 
      
     
     // printf("\n");
      //p^T * q  
      res = ddot_ac(p_ac, q_ac);
      
   
      alpha = divide(rho, res);
      
      //iter_x_ac = iter_x_ac + alpha*p
      vector_add_ac(0, alpha, p_ac, iter_x_ac, iter_x_ac);
     
      //r = r_ac - alpha*q
      vector_add_ac(0, (-alpha), q_ac, r_ac, r_ac);
      
      //re-compute the 2nd norm of r vector
      r_norm = norm_2_ac(r_ac);
      //printf("edw\n");
    }
  }
  }
  /*
  else if(cs_flag == 1){
    
      //Bi-Cg method
    if(flag == 0){
      //r_ac = b - A*iter_x_ac
      dgemv_cs_ac(0, (-1.0), C_comp, iter_x_ac, r_ac); 
      vector_add_ac(1.0, r_ac, B_ac_array_cs, r_ac);
    
      //r_tilda_ac = r_ac
      memcpy(r_tilda_ac, r_ac, size*sizeof(complex double));
      
      //norm of b vector-- only once computed
      b_norm = norm_2_ac(B_ac_array_cs);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r vector
      r_norm = norm_2_ac(r_ac);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
	
	
	preconditioner_solve_ac(M_ac, r_ac, z_ac);
	
	preconditioner_solve_ac(M_ac, r_tilda_ac, z_tilda_ac);
	
	
	//rho = r^T*z_ac
	rho = ddot_ac(z_ac, r_tilda_ac);
	
	if(fabs(rho) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
	
	if(iter == 1){
	  // p = z_ac
	  memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	  
	  // p_tilda_ac = z_tilda_ac
	  memcpy(p_tilda_ac, z_tilda_ac, (size * sizeof(complex double)));
	  
	}
	else{
	  beta = rho / rho_1;
	  
	  //p=z_ac+beta*p
	  vector_add_ac(beta, p_ac, z_ac, p_ac);
	  
	  //p_tilda_ac=z_tilda_ac+beta*p_tilda_ac
	  vector_add_ac(beta, p_tilda_ac, z_tilda_ac, p_tilda_ac);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_cs_ac(0, 1.0, C_comp, p_ac, q_ac); 
	
	//q_tilda_ac = A^T*p_tilda_ac
	dgemv_cs_ac(1, 1.0, C_comp, p_tilda_ac, q_tilda_ac); 
	
	//p^T * q  
	res = ddot_ac(p_tilda_ac, q_ac);
	
	if(fabs(res) < EPS){
	  printf("Failure with the algorithm Bi-CG\n\n");
	  
	  exit(EXIT_FAILURE);
	}
      
	alpha = rho / res;
	
	//iter_x_ac = iter_x_ac + alpha*p
	
	vector_add_ac(alpha, p_ac, iter_x_ac, iter_x_ac);
	
	
	//r = r - alpha*q
	vector_add_ac((-alpha), q_ac, r_ac, r_ac);
	
	//r_tilda_ac = r_tilda_ac - alpha*q_tilda_ac
	vector_add_ac((-alpha), q_tilda_ac, r_tilda_ac, r_tilda_ac);
	
	//re-compute the 2nd norm of r vector
	r_norm = norm_2_ac(r_ac);
	
      }
    }
    //Cg method
    else if(flag == 1){
      //r = b - A*iter_x_ac
      dgemv_cs_ac(0, (-1.0), C_comp, iter_x_ac, r_ac); 
      vector_add_ac(1.0, r_ac, B_ac_array_cs, r_ac);
    
      
      //norm of b vector-- only once computed
      b_norm = norm_2_ac(B_ac_array_cs);
      
      
      //non zero division condition
      if(b_norm == 0){
	b_norm = 1;
      }
      
      //norm of r vector
      r_norm = norm_2_ac(r_ac);

      while(((r_norm / b_norm) > itol) && (iter < MAX_ITERATION)){
	iter++;
	
	
	preconditioner_solve_ac(M_ac, r_ac, z_ac);
	
	
	//rho = r^T*z_ac
	rho = ddot_ac(r_ac, z_ac);
	
	if(iter == 1){
	  // p = z_ac
	  memcpy(p_ac, z_ac, (size * sizeof(complex double)));
	  
      
	}
	else{
	  beta = rho / rho_1;
	  
	  //p=z_ac+beta*p
	  vector_add_ac(beta, p_ac, z_ac, p_ac);
	  
	}
	
	rho_1 = rho;
	
	//q = A*p
	dgemv_cs_ac(0, 1.0, C_comp, p_ac, q_ac); 
	
	//p^T * q  
	res = ddot_ac(p_ac, q_ac);
	
      
	alpha = rho / res;
	
	//iter_x_ac = iter_x_ac + alpha*p
	
	vector_add_ac(alpha, p_ac, iter_x_ac, iter_x_ac);
	
	
	//r = r - alpha*q
	vector_add_ac((-alpha), q_ac, r_ac, r_ac);
	
	//re-compute the 2nd norm of r vector
	r_norm = norm_2_ac(r_ac);
	
      }
    }
  }*/
  
  
}

void dgemv_ac(int flag, double a, complex double *array, complex double *xx, complex double *yy){
    double value_r, value_i;
    int i,j;
    if(flag == 0){	//no traspose matrix
      for(i=0; i<size; i++){
	yy[i] = 0;
	value_r = 0;
	value_i = 0;
	for(j=0; j<size; j++){
	  value_r += ((a * creal(array[(i*size) + j])) * creal(xx[j]) - (a * cimag(array[(i*size) + j])) * cimag(xx[j]));
	  value_i += ((a * creal(array[(i*size) + j])) * cimag(xx[j]) + (a * cimag(array[(i*size) + j])) * creal(xx[j]));
	  
	}
	yy[i] = value_r + I*value_i;
      }
    }
    else if(flag == 1){	//transpose matrix
      for(i=0; i<size; i++){
	yy[i] = 0;
	value_r = 0;
	value_i = 0;
	for(j=0; j<size; j++){
	  value_r += ((a * creal(array[(j*size) + i])) * creal(xx[j]) - (a * -cimag(array[(j*size) + i])) * cimag(xx[j]));
	  value_i += ((a * creal(array[(j*size) + i])) * cimag(xx[j]) + (a * -cimag(array[(j*size) + i])) * creal(xx[j]));
	  
	}
	yy[i] = value_r + I*value_i;
      }
    }
}

void dgemv_cs_ac(int flag, double a, cs_ci *array, complex double *xx, complex double *yy){
  int j, pp;
  double value_r, value_i;
  
  for(j = 0; j<size; j++){
    yy[j] = 0;
  }
  
  if(flag == 0){	//no trapsose matrix
    for (j = 0 ; j < size ; j++){
      for (pp = array->p[j] ; pp < array->p[j+1] ; pp++){
	value_r = ((a * creal(array->x[pp])) * creal(xx[j]) - (a * cimag(array->x[pp])) * cimag(xx[j]));
	value_i = ((a * creal(array->x[pp])) * cimag(xx[j]) + (a * cimag(array->x[pp])) * creal(xx[j]));
	yy[array->i[pp]] +=  value_r + I*value_i;
      }
    }
  }
  else if(flag == 1){
    for(j = 0; j<size; j++){
      for(pp = array->p[j] ; pp < array->p[j+1] ; pp++){
	value_r = ((a * creal(array->x[pp])) * creal(xx[array->i[pp]]) - (a * -cimag(array->x[pp])) * cimag(xx[array->i[pp]]));
	value_i = ((a * creal(array->x[pp])) * cimag(xx[array->i[pp]]) + (a * -cimag(array->x[pp])) * creal(xx[array->i[pp]]));
	yy[j] += value_r + I*value_i;
      }
      
    }
  }
}


void vector_add_ac(int flag, complex double a, complex double *xx, complex double *yy, complex double *zz){
  int i;
  //double value_r, value_i;
  if(flag == 0){
    for(i=0; i<size; i++){
      zz[i] = (((creal(a) * creal(xx[i])) - (cimag(a) * cimag(xx[i]))) + creal(yy[i])) 
      + I * (((creal(a) * cimag(xx[i])) + (cimag(a) * creal(xx[i]))) + cimag(yy[i]));
    }
  }
  if(flag == 1){
    for(i=0; i<size; i++){
      zz[i] = ((a * creal(xx[i])) + creal(yy[i])) + I * ((-a * cimag(xx[i])) + cimag(yy[i]));
    }
  }
}

double norm_2_ac(complex double *xx){
  double res=0;
  int i;
  
  for(i=0; i<size; i++){
    res += (pow(creal(xx[i]), 2) + pow(cimag(xx[i]), 2)) ;
  }
  return sqrt(res);
}

complex double ddot_ac(complex double *xx, complex double *yy){
  double re = 0, im = 0;
  complex double res;
  int i;
  
  for(i=0; i<size; i++){
    
    re += ((creal(xx[i]) * creal(yy[i])) - (-cimag(xx[i]) * cimag(yy[i])));
    im += ((creal(xx[i]) * cimag(yy[i])) + (-cimag(xx[i]) * creal(yy[i])));
   
  }
  
  res = re + im * I;
  
  return res;
}

void preconditioner_solve_ac(int flag, complex double *xx, complex double *yy, complex double *zz){
  int i;
  //double value;

  if(flag == 0){
    for(i=0; i<size; i++){
	zz[i] =  divide(yy[i], xx[i]);
    }
  }
  if(flag == 1){
    for(i=0; i<size; i++){
	zz[i] =  divide(yy[i], conj(xx[i]));
      
    }
  }
}

complex double divide(complex double a, complex double b){
  complex double res;
  double value = 0;
  
  value = (creal(b) * creal(b)) + (cimag(b) * cimag(b));	//denominator
  
  res = (((creal(a) * creal(b)) + (cimag(a) * cimag(b))) / value) + 
	I * (((cimag(a) * creal(b)) - (cimag(b) * creal(a))) / value);
	
  return res;
}

/*
void double multiply(complex double a, complex double *b, complex double *y){
  int i;
  
  for(i=0; i < size; ++o){
    y[i] = ((creal(a) * creal(b[i])) - (cimag(a) * cimag(b[i]))) + I * ((creal(a) * cimag(b[i])) + (cimag(a) * creal(b[i])));
  }
}
*/

void init_cg_ac(int cs_flag, int flag){
  int i,j;
  
  //allocate memory for preconditioner M vector
  M_ac = (complex double*)calloc(size, sizeof(complex double));
  
  if(M_ac == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  if(cs_flag == 0){
  //M has the diagonal elements of MNA array
    for(i=0; i<size; i++){
      if((creal(A_ac_array[(i*size) + i]) == 0) && (cimag(A_ac_array[(i*size) + i]) == 0)){
	M_ac[i] = 1;
      }
      else{// || (cimag(A_ac_array[(i*size) + i]) != 0)){
	M_ac[i] = creal(A_ac_array[(i*size) + i]) + I * cimag(A_ac_array[(i*size) + i]);
	
      }
      
      /*\
      else{
	if(creal(A_ac_array[(i*size) + i]) == 0){
	  M_ac[i] = 1 + cimag(A_ac_array[(i*size) + i]);
	}
	else if(cimag(A_ac_array[(i*size) + i]) == 0){
	  M_ac[i] = creal(A_ac_array[(i*size) + i]) + I;
	}
      }*/
    }
  }
  else if(cs_flag == 1){
    int pp;
    
    for(j = 0; j < size; ++j){
      for(pp = C_comp->p[j]; pp < C_comp->p[j+1]; pp++){
	if(C_comp->i[pp] == j){
	  M_ac[j] = C_comp->x[pp];
	  break;
	}
      }	
      if(C_comp->i[pp] != j){
	M_ac[j] = 1; 
      }
    }	
  }
  
  
  //allocate memory for preconditioner r vector
  r_ac = (complex double*)calloc(size, sizeof(complex double));
  
  //allocate memory for preconditioner z_ac vector
  z_ac = (complex double*)calloc(size, sizeof(complex double));
  
  
  //allocate memory for preconditioner p vector
  p_ac = (complex double*)calloc(size, sizeof(complex double));
  
  
  //allocate memory for preconditioner q vector
  q_ac = (complex double*)calloc(size, sizeof(complex double));
  
  //allocate memory for preconditioner iter_x_ac vector
  iter_x_ac = (complex double*)calloc(size, sizeof(complex double));
  
  if((iter_x_ac == NULL) || (r_ac == NULL) || (z_ac == NULL) || (p_ac == NULL) || (q_ac == NULL)){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  if (flag == 0){
    //allocate memory for preconditioner r_tilda_ac vector
    r_tilda_ac = (complex double*)calloc(size, sizeof(complex double));
    
    //allocate memory for preconditioner z_tilda_ac vector
    z_tilda_ac = (complex double*)calloc(size, sizeof(complex double));
    
    
    //allocate memory for preconditioner p_tilda_ac vector
    p_tilda_ac = (complex double*)calloc(size, sizeof(complex double));
    
    
    //allocate memory for preconditioner q_tilda_ac vector
    q_tilda_ac = (complex double*)calloc(size, sizeof(complex double));
    
    if((r_tilda_ac == NULL) || (z_tilda_ac == NULL) || (p_tilda_ac == NULL) || (q_tilda_ac == NULL)){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
    }
  }
  
}

void printIter_x_ac(int flag){
  int i;
  
  if(flag == 0){
    
    printf("\n*****Bi-CG METHOD*****\n\n");
    printf("The Bi-CG algorith ran %d times until solution...\n\n", iter);
    //print the x vector
    printf ("x = \n\n");
    
    for(i=0; i<size; i++){
      printf("%e + j*(%e)\n", creal(iter_x_ac[i]), cimag(iter_x_ac[i]));
    }
    
    printf("\n\n");
  }
  else if(flag == 1){
    printf("\n*****CG METHOD*****\n\n");
    printf("The CG algorith ran %d times until solution...\n\n", iter);
    //print the x vector
    printf ("x = \n\n");
    
    for(i=0; i<size; i++){
      printf("%e + j*(%e)\n", creal(iter_x_ac[i]), cimag(iter_x_ac[i]));
    }
    
    printf("\n\n");
  }
}


void freeIter_ac(){
  
  free(iter_x_ac);
  
  free(r_ac);
  
  free(z_ac);
  
  free(q_ac);
  
  free(p_ac);
  
  free(M_ac);
  
  free(r_tilda_ac);
  
  free(z_tilda_ac);
  
  free(q_tilda_ac);
  
  free(p_tilda_ac);
  
    
}