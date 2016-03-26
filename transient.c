#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "hashtable.h"
#include "mna.h"
#include "mna_sparse.h"
#include "parser.h"
#include "factorization.h"
#include "iter_method.h"
#include "csparse.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "transient.h"

#define PI 3.14159265

int size, k;

double *A_temp;
double *x_prev;
double *b_temp;


double *e_current = NULL;
double *e_prev = NULL;

gsl_vector *e_g;
gsl_vector *e_prev_g;


void transient(int n, int m2, double time_step, double finish_time, int m_flag, int f_flag, int i_flag, int s_flag, double itol){
  //double k;
  double t;
  size = n + m2;
  
  if(s_flag == 0){
    if(i_flag == 0){
      if(m_flag == 0){
	trapezoidal(time_step, finish_time, f_flag);
      }
      else if(m_flag == 1){
	backward_euler(time_step, finish_time, f_flag);
      }
    }
    else if(i_flag == 1){
      if(m_flag == 0){
	trapezoidal_iter(n, m2, time_step, finish_time, f_flag, s_flag, itol);
      }
      else if(m_flag == 1){
	backward_euler_iter(n, m2, time_step, finish_time, f_flag, s_flag, itol);
      }
    }
  }
  else if(s_flag == 1){
    if(i_flag == 0){
      if(m_flag == 0){
	trapezoidal_cs(time_step, finish_time, f_flag);
      }
      else if(m_flag == 1){
	backward_euler_cs(time_step, finish_time, f_flag);
      }
    }
    else if(i_flag == 1){
      if(m_flag == 0){
	trapezoidal_iter_cs(n, m2, time_step, finish_time, f_flag, s_flag, itol);
      }
      else if(m_flag == 1){
	backward_euler_iter_cs(n, m2, time_step, finish_time, f_flag, s_flag, itol);
	
      }
    }
  }
  
  
  //tran_printf();
  
  
}

void backward_euler(double time_step, double finish_time, int flag){
    int i, m, j=0, l;
    
    printf("*****Backward Euler method*****\n\n");
    
    k = finish_time / time_step;
    
    /*
    double nn = fmod(finish_time, time_step);
    printf("%g\n", nn);
    
    if( fmod(finish_time, time_step) != 0){
      printf("okk\n");
      k++;
    }
    */
    
    A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    //G + (1/h * C)
    for(i=0; i < size*size; i++){
      A_temp[i] = A_array[i] + (1/time_step) * C_array[i];
    }
    
    
    //copy all elements of A temp array to temp array
    for(i=0; i < size*size; i++){
      tmp_array[i] = A_temp[i];
    }
      
      
    if(flag == 0){
      lu_Decomp();
  
    }
    else if(flag == 1){  
      chol_Decomp();
    }
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = gsl_vector_get(x, m);
      j++;
    }
    
    for(l = 1; l <= k; l++){
    
      
	compute_E(l*time_step);
      
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = gsl_vector_get(x, i);
      }
      
      //e(t_k) + 1/h * C * x(t_k - 1)
      dgemv_tran(0, 1/time_step, C_array, x_prev, b_temp);
      vector_add_tran(1.0, b_temp, e_current, b_temp);
      
    
      //copy all elements of b temp vector to temp2 array
      for(i=0; i < size; i++){
	tmp2_array[i] = b_temp[i];
      }
      
      if(flag == 0){
	lu_solve();
      }
      else if(flag == 1){
	chol_solve();
      }
      
      //keep the rest n elements of x vector
	for(m=0; m < size; m++){
	  tran_array[j] = gsl_vector_get(x, m);
	  j++;
	}
	
     
  }
}

void backward_euler_cs(double time_step, double finish_time, int flag){
    int i, m, j=0, l;
    
    printf("*****Backward Euler method for sparse matrices*****\n\n");
    
    k = finish_time / time_step;
    
    /*
    if( fmod(finish_time, time_step) != 0){
      printf("ok\n");
      //k++;
    }
    */
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    //b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    

    for(i=0; i<size; i++){
      tran_array[j] = B_array_cs_temp[i];
      j++;
    }
    
    C_comp = cs_compress(Asparse);	//create compressed column matrix
    
    cs_dupl(C_comp);
    
    C_tilda_comp = cs_compress(Csparse);	//create compressed column matrix
    
    cs_dupl(C_tilda_comp);
  
     //G + (1/h * C)
    C_comp = cs_add(C_comp, C_tilda_comp, 1, 1/time_step);

    factorization_cs(flag);   
    
    for(l = 1; l <= k; l++){
    
      compute_E(l*time_step);
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = B_array_cs_temp[i];
      }
      

      //e(t_k) + 1/h * C * x(t_k - 1)
      dgemv_tran_cs(0, 1/time_step, C_tilda_comp, x_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_current, B_array_cs);
      
   

      if(flag == 0){
	lu_solve_cs(1);
      }
      else if(flag == 1){
	chol_solve_cs(1);
      }
      
                 

      //keep the rest n elements of x vector
	for(m=0; m < size; m++){
	  tran_array[j] =B_array_cs_temp[m];
	  j++;
	}
	
     
  }
  
 
    
    
}


void backward_euler_iter(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol){
    int i, m, j=0, l;
   
    
    printf("*****Backward Euler iterative method*****\n\n");
    
    k = finish_time / time_step;
    
    /*
    if( fmod(finish_time, time_step) != 0){
      k++;
    }
    */
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    //b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    
    //G + (1/h * C)
    for(i=0; i < size*size; i++){
      A_array[i] += ((1/time_step) * C_array[i]);
    }
    
    //M has the diagonal elements of new A_array
    for(i=0; i<size; i++){
      if(A_array[(i*size) + i] == 0){
	M[i] = 1;
      }
      else{
	M[i] =  A_array[(i*size) + i];
      }
    }
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = iter_x[m];
      j++;
    }
    
    for(l = 1; l <= k; l++){
    
      compute_E(l*time_step);
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] =  iter_x[i];
      }
      
      //e(t_k) + 1/h * C * x(t_k - 1)
      dgemv_tran(0, 1/time_step, C_array, x_prev, B_array);
      vector_add_tran(1.0, B_array, e_current, B_array);
      
    
      
      if(flag == 0){
	cg_method(cs_flag, 0, 1, itol, n, m2);
      }
      else if(flag == 1){
	cg_method(cs_flag, 1, 1, itol, n, m2);
      }
      
      //keep the rest n elements of x vector
	for(m=0; m < size; m++){
	  tran_array[j] = iter_x[m];
	  j++;
	}
	
     
    }
}


void backward_euler_iter_cs(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol){
    int i, m, j=0, l;
   
    
    printf("*****Backward Euler iterative method for sparse matrices*****\n\n");
     
    k = finish_time / time_step;
    
    /*
    if( fmod(finish_time, time_step) != 0){
      k++;
    }
    */
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    //b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
   
    C_comp = cs_compress(Asparse);	//create compressed column matrix
    
    cs_dupl(C_comp);
    
    C_tilda_comp = cs_compress(Csparse);	//create compressed column matrix
    
    cs_dupl(C_tilda_comp);
  
     //G + (1/h * C)
    C_comp = cs_add(C_comp, C_tilda_comp, 1, 1/time_step);

    
    //M has the diagonal elements of new A_array
     int pp;
    
    for(i = 0; i < size; ++i){
      for(pp = C_comp->p[i]; pp < C_comp->p[i+1]; pp++){
	if(C_comp->i[pp] == i){
	  M[i] = C_comp->x[pp];
	  break;
	}
      }	
      if(C_comp->i[pp] != i){
	M[i] = 1; 
      }
    }	
  
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = iter_x[m];
      j++;
    }
    
    for(l = 1; l <= k; l++){
    
      compute_E(l*time_step);
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] =  iter_x[i];
      }
      
      //e(t_k) + 1/h * C * x(t_k - 1)
      dgemv_tran_cs(0, 1/time_step, C_tilda_comp, x_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_current, B_array_cs);
      
    
      
      if(flag == 0){
	cg_method(cs_flag, 0, 1, itol, n, m2);
      }
      else if(flag == 1){
	cg_method(cs_flag, 1, 1, itol, n, m2);
      }
      
      //keep the rest n elements of x vector
	for(m=0; m < size; m++){
	  tran_array[j] = iter_x[m];
	  j++;
	}
	
     
    }
}


void trapezoidal(double time_step, double finish_time, int flag){
    int i, m, j=0, l;
    
    printf("*****Trapezoidal method*****\n\n");
   
    k = finish_time / time_step;
    /*
    if( fmod(finish_time, time_step) != 0){
      k++;
    }
    */
    A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    e_prev = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    
    //G + (2/h * C)
    for(i=0; i < size*size; i++){
      A_temp[i] = A_array[i] + (2/time_step) * C_array[i];
    }
    
    //copy all elements of A temp array to temp array
    for(i=0; i < size*size; i++){
      tmp_array[i] = A_temp[i];
    }
    
    //G - (2/h * C)
    for(i=0; i < size*size; i++){
      C_array[i] = A_array[i] - (2/time_step) * C_array[i];
    }
      
    if(flag == 0){
      lu_Decomp();
  
    }
    else if(flag == 1){  
      chol_Decomp();
    }
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = gsl_vector_get(x, m);
      j++;
    }
    
    //compute e_prev for the first time
    compute_E(0);
    
    for(m=0; m < size; m++){
      e_prev[m] = e_current[m];
    }
    
    for(l = 1; l <= k; l++){
      
      //compute e_current
      compute_E(l*time_step);
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = gsl_vector_get(x, i);
      }
      
      //e(t_k) + e(t_k - 1) - (G - 2/h * C) * x(t_k - 1)
      dgemv_tran(0, (-1.0), C_array, x_prev, b_temp);
      vector_add_tran(1.0, b_temp, e_prev, b_temp);
      vector_add_tran(1.0, b_temp, e_current, b_temp);
      
    
      //copy all elements of b temp vector to temp2 array
      for(i=0; i < size; i++){
	tmp2_array[i] = b_temp[i];
      }
      
      if(flag == 0){
	lu_solve();
      }
      else if(flag == 1){
	chol_solve();
      }
      
      //keep the rest n elements of x vector
      for(m=0; m < size; m++){
	tran_array[j] = gsl_vector_get(x, m);
	j++;
      }
      
      for(m=0; m < size; m++){
	e_prev[m] = e_current[m];
      }
      
    }
  
}


void trapezoidal_cs(double time_step, double finish_time, int flag){
    int i, m, j=0, l;
    
    printf("*****Trapezoidal method for sparse matrices*****\n\n");
   
    k = finish_time / time_step;
    
    ///if( fmod(finish_time, time_step) != 0){
     // k++;
   // }
    
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
   // b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    e_prev = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    //keep the elements for x_t0
    for(i=0; i<size; i++){
      tran_array[j] = B_array_cs_temp[i];
      j++;
    }
    
    
    C_comp = cs_compress(Asparse);	//create compressed column matrix
    
    cs_dupl(C_comp);
    
    C_tilda_comp = cs_compress(Csparse);	//create compressed column matrix
    
    cs_dupl(C_tilda_comp);
    
    
    
     //G + (2/h * C)
    C_comp = cs_add(C_comp, C_tilda_comp, 1, 2/time_step);
    
    
    //G - (2/h * C)
    C_tilda_comp = cs_add(C_comp, C_tilda_comp, 1, (-4/time_step));
    
    
    factorization_cs(flag); 
    
    
    //compute e_prev for the first time
    compute_E(0);
    
    for(m=0; m < size; m++){
      e_prev[m] = e_current[m];
    }
    
    //printf("re\n");
    for(l = 1; l <= k; l++){
      
      //compute e_current
      compute_E(l*time_step);
      
      
     //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = B_array_cs_temp[i];
      }
      
      
      //e(t_k) + e(t_k - 1) - (G - 2/h * C) * x(t_k - 1)
      dgemv_tran_cs(0, (-1.0), C_tilda_comp, x_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_current, B_array_cs);
      
      
      if(flag == 0){
	lu_solve_cs(1);
      }
      else if(flag == 1){
	chol_solve_cs(1);
      }
      
      //keep the rest n elements of x vector
      for(m=0; m < size; m++){
	tran_array[j] = B_array_cs_temp[m];
	j++;
      }
    
      for(m=0; m < size; m++){
	e_prev[m] = e_current[m];
      }
      
    }
  
}


void trapezoidal_iter(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol){
    int i, m, j=0, l;
    
    printf("*****Trapezoidal iterative method*****\n\n");
   
    k = finish_time / time_step;
    
    //if( fmod(finish_time, time_step) != 0){
      //k++;
    //}
    
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    //b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    e_prev = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    
    //G + (2/h * C)
    for(i=0; i < size*size; i++){
      A_array[i] += (2/time_step) * C_array[i];
    }
    
    //M has the diagonal elements of new A_array
    for(i=0; i<size; i++){
      if(A_array[(i*size) + i] == 0){
	M[i] = 1;
      }
      else{
	M[i] =  A_array[(i*size) + i];
      }
    }
    
    //G - (2/h * C)
    for(i=0; i < size*size; i++){
      C_array[i] = A_array[i] - ((4/time_step) * C_array[i]);
    }
      
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = iter_x[m];
      j++;
    }
    
     //compute e_prev for the first time
    compute_E(0);
    
    for(m=0; m < size; m++){
      e_prev[m] = e_current[m];
    }
   
    for(l = 1; l <= k; l++){
      
      //compute e_current
      compute_E(l*time_step);
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = iter_x[i];
      }
      
      //e(t_k) + e(t_k - 1) - (G - 2/h * C) * x(t_k - 1)
      dgemv_tran(0, (-1.0), C_array, x_prev, B_array);
      vector_add_tran(1.0, B_array, e_prev, B_array);
      vector_add_tran(1.0, B_array, e_current, B_array);
      
    
      
      if(flag == 0){
	cg_method(cs_flag, 0, 1, itol, n, m2);
      }
      else if(flag == 1){
	cg_method(cs_flag, 1, 1, itol, n, m2);
      }
      
      //keep the rest n elements of x vector
      for(m=0; m < size; m++){
	tran_array[j] = iter_x[m];
	j++;
      }
    
      for(m=0; m < size; m++){
	e_prev[m] = e_current[m];
      }
      
    }
  
}


void trapezoidal_iter_cs(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol){
    int i, m, j=0, l;
    
    printf("*****Trapezoidal iterative method for sparse matrices*****\n\n");
   
    k = finish_time / time_step;
    
    //if( fmod(finish_time, time_step) != 0){
      //k++;
    //}
    
    //A_temp = (double*)calloc(size*size, sizeof(double));
    x_prev = (double*)calloc(size, sizeof(double));
    //b_temp = (double*)calloc(size, sizeof(double));
    e_current = (double*)calloc(size, sizeof(double));
    e_prev = (double*)calloc(size, sizeof(double));
    tran_array = (double*)calloc((k+1) * size , sizeof(double));
    
    
    C_comp = cs_compress(Asparse);	//create compressed column matrix
    
    cs_dupl(C_comp);
    
    C_tilda_comp = cs_compress(Csparse);	//create compressed column matrix
    
    cs_dupl(C_tilda_comp);
    
    //G + (2/h * C)
    C_comp = cs_add(C_comp, C_tilda_comp, 1, 2/time_step);
    
    
    //G - (2/h * C)
    C_tilda_comp = cs_add(C_comp, C_tilda_comp, 1, (-4/time_step));
    
    
      //M has the diagonal elements of new A_array
     int pp;
    
    for(i = 0; i < size; ++i){
      for(pp = C_comp->p[i]; pp < C_comp->p[i+1]; pp++){
	if(C_comp->i[pp] == i){
	  M[i] = C_comp->x[pp];
	  break;
	}
      }	
      if(C_comp->i[pp] != i){
	M[i] = 1; 
      }
    }	
  
      
    
    //keep the first n elements of x vector for x0
    for(m=0; m < size; m++){
      tran_array[j] = iter_x[m];
      j++;
    }
    
     //compute e_prev for the first time
    compute_E(0);
    
    for(m=0; m < size; m++){
      e_prev[m] = e_current[m];
    }
   
    for(l = 1; l <= k; l++){
      
      //compute e_current
      compute_E(l*time_step);
      
      
      //take the x_prev from the x solution factorization
      for(i=0; i < size; i++){
	x_prev[i] = iter_x[i];
      }
      
      //e(t_k) + e(t_k - 1) - (G - 2/h * C) * x(t_k - 1)
      dgemv_tran_cs(0, (-1.0), C_tilda_comp, x_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_prev, B_array_cs);
      vector_add_tran(1.0, B_array_cs, e_current, B_array_cs);
      
    
      
      if(flag == 0){
	cg_method(cs_flag, 0, 1, itol, n, m2);
      }
      else if(flag == 1){
	cg_method(cs_flag, 1, 1, itol, n, m2);
      }
      
      
      //keep the rest n elements of x vector
      for(m=0; m < size; m++){
	tran_array[j] = iter_x[m];
	j++;
      }
    
      for(m=0; m < size; m++){
	e_prev[m] = e_current[m];
      }
      
    }
  
}


void compute_E(double t){
  element *tmp;
  pwl *p_tmp, *p_prev;
  double value = 0;
  int i,  pos, neg, l = 0;
  
  tmp = head;
  
  memset(e_current, 0, size*sizeof(double));
  
  while(tmp!=NULL){
    if(tmp->type_id== 'v'){
      Volt_Source* temp = (Volt_Source*)(tmp->type);
      
      if(temp->t_flag == 1){
	
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  
	  if(t <= t_temp->t1){
	    value = t_temp->i1;
	  }
	  else{
	     //if t is between struct's t1 kai head's t1
	      if(t <= p_head->t1){
		value = linear_value(t_temp->t1, p_head->t1, t_temp->i1, p_head->i1, t);
		tmp = tmp->next;
		e_current[n+l] += value;
      
		l++;
		continue;
	      }
	  
	    p_tmp = p_head;
	    p_prev = p_tmp;
	    p_tmp = p_tmp->next;
	    
	    while(p_tmp!=NULL){
	     
	      if(t <= p_tmp->t1){
		value = linear_value(p_prev->t1, p_tmp->t1, p_prev->i1, p_tmp->i1, t);
		break;		
	      }
	      p_prev = p_tmp;
	      p_tmp = p_tmp->next;
	    }
	    if(p_tmp == NULL){
	      value = p_prev->i1;
	    }
	    
	  }
	}
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  if(t <= t_temp->td1){
	    value = t_temp->i1;
	  }
	  else if(t <= t_temp->td2){
	    value = t_temp->i1 + (t_temp->i2 - t_temp->i1)*(1 - exp(-(t - t_temp->td1) / t_temp->tc1));
	  }
	  else{
	    value = t_temp->i1 + (t_temp->i2 - t_temp->i1) * (exp(-(t - t_temp->td2) / t_temp->tc2) - exp(-(t - t_temp->td1) / t_temp->tc1));
	  }
	}
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  if(t <= t_temp->td){
	    value = t_temp->i1 + t_temp->ia * sin(2 * PI * t_temp->ph / 360); 
	  }
	  else{
	    value = t_temp->i1 + t_temp->ia * sin(2 * PI * t_temp->fr * (t - t_temp->td) + (2 * PI * t_temp->ph / 360)) * exp(-(t - t_temp->td) * t_temp->df);
	  }
	}
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  
	  int kk = (t - t_temp->td)/ t_temp->per;
	  if(kk<0){
	    kk = 0;
	  }
	  
	  t -= kk * t_temp->per;
	  
	  if(t <= t_temp->td){
	    value = t_temp->i1;
	  }
	  else if(t <= t_temp->td + t_temp->tr){
	    value = linear_value(t_temp->td, t_temp->td+t_temp->tr, t_temp->i1, t_temp->i2, t);
	  }
	  else if(t <= t_temp->td + t_temp->tr + t_temp->pw){
	    value = t_temp->i2;
	  }
	  else if(t <= t_temp->td + t_temp->tr + t_temp->pw + t_temp->tf){
	    value = linear_value(t_temp->td + t_temp->tr + t_temp->pw, t_temp->td + t_temp->tr + t_temp->pw + t_temp->tf , t_temp->i2, t_temp->i1, t );
	  }
	  else{
	    value = t_temp->i1;
	  }
	}
      }
      else{
	value = temp->value;
      }
      
      
      e_current[n+l] +=  value;
      
      l++;
      
    }
    else if(tmp->type_id== 'i'){
      Power_Source* temp = (Power_Source*)(tmp->type);
      
      if(temp->t_flag == 1){
	
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  
	  if(t <= t_temp->t1){
	    value = t_temp->i1;
	  }
	  else{
	     //if t is between struct's t1 kai head's t1
	      if(t <= p_head->t1){
		value = linear_value(t_temp->t1, p_head->t1, t_temp->i1, p_head->i1, t);
		tmp=tmp->next;
		
		//if positive node == 0 do nothing
		if(strcmp(temp->pos_node, "0") != 0){
		  pos = search_element_id(temp->pos_node);
		  
		  e_current[pos] -= value;
		  
		}
		//if negative node == 0 do nothing
		if(strcmp(temp->neg_node, "0") != 0){
		  neg = search_element_id(temp->neg_node);
		  
		  e_current[neg] +=  value;
		}
      
		continue;
	      }
	    
	    p_tmp = p_head;
	    p_prev = p_tmp;
	    p_tmp = p_tmp->next;
	    
	    while(p_tmp!=NULL){
	     
	      if(t <= p_tmp->t1){
		value = linear_value(p_prev->t1, p_tmp->t1, p_prev->i1, p_tmp->i1, t);
		break;		
	      }
	      p_prev = p_tmp;
	      p_tmp = p_tmp->next;
	     	      
	    }
	    if(p_tmp == NULL){
	      value = p_prev->i1;
	    }
	  }
	}
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  if(t <= t_temp->td1){
	    value = t_temp->i1;
	  }
	  else if(t <= t_temp->td2){
	    value = t_temp->i1 + (t_temp->i2 - t_temp->i1)*(1 - exp(-(t - t_temp->td1) / t_temp->tc1));
	  }
	  else{
	    value = t_temp->i1 + (t_temp->i2 - t_temp->i1) * (exp(-(t - t_temp->td2) / t_temp->tc2) - exp(-(t - t_temp->td1) / t_temp->tc1));
	  }
	}
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  if(t <= t_temp->td){
	    value = t_temp->i1 + t_temp->ia * sin(2 * PI * (t_temp->ph / 360)); 
	  }
	  else{
	    value = t_temp->i1 + t_temp->ia * sin(2 * PI * t_temp->fr * (t - t_temp->td) + (2 * 3.14 * t_temp->ph / 360)) * exp(-(t - t_temp->td) * t_temp->df);
	  }
	}
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  
	  int kk = (t - t_temp->td)/ t_temp->per;
	  if(kk<0){
	    kk = 0;
	  }
	  
	  t -= kk * t_temp->per;
	  
	  
	  if(t <= t_temp->td){
	    value = t_temp->i1;
	  }
	  else if(t <= t_temp->td + t_temp->tr){
	    value = linear_value(t_temp->td  , t_temp->td+t_temp->tr, t_temp->i1, t_temp->i2, t);
	  }
	  else if(t <= t_temp->td + t_temp->tr + t_temp->pw){
	    value = t_temp->i2;
	  }
	  else if(t <= t_temp->td + t_temp->tr + t_temp->pw + t_temp->tf){
	    value = linear_value(t_temp->td + t_temp->tr + t_temp->pw, t_temp->td + t_temp->tr + t_temp->pw + t_temp->tf, t_temp->i2, t_temp->i1, t);
	  }
	  else{
	    value = t_temp->i1;
	  }
	}
      }
      else{
	value = temp->value;
      }
      
      //if positive node == 0 do nothing
      if(strcmp(temp->pos_node, "0") != 0){
	pos = search_element_id(temp->pos_node);
	
	e_current[pos] -= value;
	
      }
      //if negative node == 0 do nothing
      if(strcmp(temp->neg_node, "0") != 0){
	neg = search_element_id(temp->neg_node);
	
	e_current[neg] += value;
      }
            
    }
    else if(tmp->type_id== 'l'){
      l++;
    }
    
    tmp = tmp->next;
  }
  
}


double linear_value(double x1, double x2, double y1, double y2, double x){
    double y, yy;
    
    y = (y1 - y2) / (x1 - x2);
  
    yy = y * (x - x2) + y2;
  
    return yy ;
}

void dgemv_tran(int flag, double a, double *array, double *xx, double *yy){
    
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

void dgemv_tran_cs(int flag, double a, cs *array, double *xx, double *yy){
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



void vector_add_tran(double a, double *xx, double *yy, double *zz){
    int i;
    
    for(i=0; i<size; i++){
      zz[i] = (a *xx[i]) + yy[i];
    }
}

void freeTran(){
  
    free(A_temp);
    
    free(b_temp);
    
    free(x_prev);
    
    free(e_current);
    
    free(e_prev);
    
    free(tran_array);
}

void tran_printf(){
    int i, j;
    
    
    printf("Transient analysis: \n\n");
    
    
    //print all the x vectors for the transient analysis
    for(i=0; i < (k+1); i++){
      for(j=0; j < size; j++){
	printf("x%d_%d = %.5lf\t", (i+1), (j+1), tran_array[(i*size) + j]);
      }
      printf("\n");
    }
    
    printf("\n\n");
}