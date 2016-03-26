/* Georgiadis Xristos 588
 * Stabolos Nikolaos 594
 * Vasileiadou Pelagia 1054
 * Pappas-Vafias Panagiotis 1106
*/
/* ######################## STAGE 1 ######################## */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex_math.h>
#include <complex.h>
#include "parser.h"
#include "hashtable.h"
#include "mna.h"
#include "factorization.h"
#include "dc_sweep.h"
#include "options.h"
#include "iter_method.h"
#include "iter_method_ac.h"
#include "mna_sparse.h"
#include "cs.h"
#include "transient.h"
#include "csparse.h"
#include "ac_analysis.h"
#define PI 3.14159265

element *curr = NULL;
dc *dc_tmp = NULL;
plot *plot_tmp = NULL;


int main(int argc, char *argv[]){
  FILE *fp;
  char *line = NULL;
  size_t read;
  size_t len;
  char *token;
  char *value;
  double itol, tmp, tmp1, time_step, finish_time;
  char *input_variable=NULL;
  char del[] = " \t\n\r";
  int line_num=0, ground_flag=0, m2=0, dc_sweep_flag = 0, plot_flag = 0, ac_flag = 0;
  int factor_flag , tran_flag, iter_flag, sparse_flag, method_flag, m, cholesky_flag = 0;
  int id = 0, n_z=0, C_n_z=0, j=0, i, ii, ac_n_z=0;	//counter for the numbers of .plot 
  int lin_flag = 0, log_flag =0;
  struct timespec  tv1, tv2;
   char *ac_method;
   int points;
   double start_freq, end_freq, tmp_freq = 0;
   double value_r, value_i;
   
	  
  //open netlist file
  fp = fopen(argv[1], "r");
  if(fp==NULL){
    printf("File %s not exist\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  
  printf("\n");
  
  //clock starts
  clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);
  
  while((read = getline(&line, &len, fp)) != -1){
    
    //integer for the line number
    line_num++;
   
    //get the first token of the line
    token = strtok(line,del);
    
    if(token[0] != '*'){ //check if it is a comment
      if((token[0]== 'v') || (token[0] == 'V')){
	Volt_Source *vSource;
	vSource = (Volt_Source*)malloc(sizeof(Volt_Source));
	
	//check for memory error
	if(vSource == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	volt_info(vSource, token, &ground_flag);// Volt Source
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = vSource;
	curr->type_id = 'v';
	
	//insert nodes in hashtable
	Volt_Source* temp = (Volt_Source*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
	
	//check for non_zero elements in mna
	if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	  n_z += 4;
	  ac_n_z += 4;//+2 for -jωl
	}
	else{
	  n_z += 2;
	  ac_n_z += 2;//+2 for -jωl
	}
	
	//counter for team_2 objects
	m2++;
	
	cholesky_flag = 1;	//flag for the existence of Volt Source or Inductance. 
      }
      else if((token[0]== 'i') || (token[0] == 'I')){
	Power_Source *pSource;
	pSource = (Power_Source*)malloc(sizeof(Power_Source));
	
	//check for memory error
	if(pSource == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	power_info(pSource, token, &ground_flag);//Power source
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = pSource;
	curr->type_id = 'i';
	
	//insert nodes in hashtable
	Power_Source* temp = (Power_Source*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
	
      }
      else if((token[0]== 'r') || (token[0] == 'R')){
	Resistance *resistance;
	resistance = (Resistance*)malloc(sizeof(Resistance));
	
	//check for memory error
	if(resistance == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	resistance_info(resistance, token, &ground_flag);//Resistance
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = resistance;
	curr->type_id = 'r';
	
	//insert nodes in hashtable
	Resistance* temp = (Resistance*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
	
	//check for non_zero elements in mna
	if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	  n_z += 4;
	  ac_n_z += 4;//+4 for -jωl
	}
	else{
	  n_z += 1;
	  ac_n_z += 1;//+1 for -jωl
	}
	
      }
      else if((token[0]== 'c') || (token[0] == 'C')){
	Capacity *capacity;
	capacity = (Capacity*)malloc(sizeof(Capacity));
	
	//check for memory error
	if(capacity == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	capacity_info(capacity, token, &ground_flag);//Capacity
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = capacity;
	curr->type_id = 'c';
	
	//insert nodes in hashtable
	Capacity* temp = (Capacity*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
	
	//check for non_zero elements in mna
	if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	  C_n_z += 4;
	  ac_n_z += 4;//+4 for -jωl
	}
	else{
	  C_n_z += 1;
	  ac_n_z += 1;//+1 for -jωl
	}
	
      }
      else if((token[0]== 'l') || (token[0] == 'L')){
	Inductance *inductance;
	inductance = (Inductance*)malloc(sizeof(Inductance));
	
	//check for memory error
	if(inductance == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	inductance_info(inductance, token, &ground_flag);//Inductance
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = inductance;
	curr->type_id = 'l';
	
	//insert nodes in hashtable
	Inductance* temp = (Inductance*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
	
	//check for non_zero elements in mna
	if((strcmp(temp->pos_node, "0") != 0) && (strcmp(temp->neg_node, "0") != 0)){
	  n_z += 4;
	  ac_n_z += 4;
	}
	else{
	  n_z += 2;
	  ac_n_z += 2;
	}
	
	ac_n_z++; //+1 for -jωl
	
	//counter for non_zero elements in C tilda matrix
	C_n_z++;
	
	//counter for team_2 objects
	m2++;
	
	cholesky_flag = 1;	//flag for the existence of Volt Source or Inductance. 
      }
      else if((token[0]== 'd') || (token[0] == 'D')){
	Diode *diode;
	diode = (Diode*)malloc(sizeof(Diode));
	
	//check for memory error
	if(diode == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	diode_info(diode, token, &ground_flag);
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = diode;
	curr->type_id = 'd';
	
	//insert nodes in hashtable
	Diode* temp = (Diode*)(curr->type);
	insert(temp->pos_node, token[0]);
	insert(temp->neg_node, token[0]);
      }
      else if((token[0]== 'm') || (token[0] == 'M')){
	MOS_transistor *mos_tran;
	mos_tran = (MOS_transistor*)malloc(sizeof(MOS_transistor));
	
	//check for memory error
	if(mos_tran == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	mos_info(mos_tran, token, &ground_flag);
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = mos_tran;
	curr->type_id = 'm';
	
	//insert nodes in hashtable
	MOS_transistor* temp = (MOS_transistor*)(curr->type);
	insert(temp->drain, token[0]);
	insert(temp->gate, token[0]);
	insert(temp->source, token[0]);
	insert(temp->body, token[0]);
      }
      else if((token[0]== 'q') || (token[0] == 'Q')){
	BJT_transistor *bjt_tran;
	bjt_tran = (BJT_transistor*)malloc(sizeof(BJT_transistor));
	
	//check for memory error
	if(bjt_tran == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	bjt_info(bjt_tran, token, &ground_flag);
	
	//insert to the list
	curr = insert_element(curr);
	curr->type = bjt_tran;
	curr->type_id = 'q';
	
	//insert nodes in hashtable
	BJT_transistor* temp = (BJT_transistor*)(curr->type);
	insert(temp->collector, token[0]);
	insert(temp->base, token[0]);
	insert(temp->emitter, token[0]);
      }
      else if(token[0] == '.'){
	if((strcmp(&token[1], "DC") == 0) || (strcmp(&token[1], "dc") == 0)){
	  token = strtok(NULL, del);
	  
	  dc_sweep_flag = dc_option(token); //dc sweep for Volt or Power
	  
	}	
	else if((strcmp(&token[1], "TRAN") == 0) || (strcmp(&token[1], "tran") == 0)){
	  
	  //transient analysis
	  token = strtok(NULL, del);
	  
	  time_step = strtod(token, &value);
	  
	  token = strtok(NULL, del);
	  
	  finish_time = strtod(token, &value);
	  
	  read = getline(&line, &len, fp);
	  
	  //get the first token of the line. Plot must be after tran in this case
	  token = strtok(line,del);
	  
	  if(token[0] != '.'){
	    printf("ERROR after .TRAN...\n");
	    
	    
	    fclose(fp);
	    if(line){
	      free(line);
	    }
	
	    //free list
	    free_list();

	    //free hashtable
	    freeHashTable();
	    exit(EXIT_FAILURE);
	  }
	  
	  tran_plot(token);
	  
	  tran_flag = 1;
	}
	else if((strcmp(&token[1], "AC") == 0) || (strcmp(&token[1], "ac") == 0)){
	 //ac analysis
	  
	  token = strtok(NULL, del);
	  
	  //LIN or LOG
	  ac_method = (char*)malloc((strlen(token)+1)*sizeof(char));
	  
	  
	  if((strncmp(token, "LIN", 3) == 0) || (strncmp(token, "lin", 3) == 0)){
	    lin_flag = 1;
	  }
	  else if((strncmp(token, "LOG", 3) == 0) || (strncmp(token, "log", 3) == 0)){
	    log_flag = 1;
	  }
	  
	  strcpy(ac_method, token);
	  
	  token = strtok(NULL, del);
	  
	  //points
	  points = atoi(token);
	  
	  token = strtok(NULL, del);
	  
	  //start_freq
	  start_freq = strtod(token, &value);
	  
	  token = strtok(NULL, del);
	  
	  //end_freq
	  end_freq = strtod(token, &value);
	  
	  read = getline(&line, &len, fp);
	  
	  //get the first token of the line. plot must be under ac in this case
	  token = strtok(line,del);
	  
	  if(token[0] != '.'){
	    printf("ERROR after .AC...\n");
	    
	    fclose(fp);
	    if(line){
	      free(line);
	    }
	
	    //free list
	    free_list();

	    //free hashtable
	    freeHashTable();
	    exit(EXIT_FAILURE);
	  }
	  
	  ac_plot(token);
	  
	  ac_flag = 1;
	  
	}
	else if((strcmp(&token[1], "PLOT") == 0) || (strcmp(&token[1], "plot") == 0)){
	  token = strtok(NULL, del);
	  
	  plot_option(token, id);
	  
	  id++;		//we assume that we have the same number of .DC and .PLOT in the netlist
	  plot_flag = 1;
	}
	else if((strcmp(&token[1], "OPTIONS") == 0) || (strcmp(&token[1], "options") == 0)){
	  all_options(token);
	  
	  factor_flag = factor_f;
	  iter_flag = iter_f;
	  itol = itol_t;
	  sparse_flag = sparse_f;
	  method_flag = method_f;
	}
      }
      else{
	printf("ERROR in line %d of the input file: '%s' element cannot be recognized !!!\n", line_num, token);
	
	fclose(fp);
	if(line){
	  free(line);
	}
	
	//free list
	free_list();

	//free hashtable
	freeHashTable();
	
	exit(EXIT_FAILURE);
      }
    }
  }
  
  /* This is the end of the main computation. Take the end time,  *
	 * calculate the duration of the computation and report it. 	*/
    clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);
  
    printf ("\n\nTotal time for reading %d lines in netlist = %10g seconds\n\n", line_num, 
		  (double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
		  (double) (tv2.tv_sec - tv1.tv_sec));
  
    
  fclose(fp);
  if(line)
    free(line);
  
  
  if(ground_flag != 1) {
    printf("The circuit cannot be simulated. Ground node does NOT exist!\n\n");
    
    fclose(fp);
    if(line){
      free(line);
    }
    
    //free list
    free_list();

    //free hashtable
    freeHashTable();
    
    freeList();

    freeList2();
  
    freeList3();
    
    exit(EXIT_FAILURE);
  }

  /* ######################## PRINT ######################## */
  if(sparse_flag == 0){
    if(iter_flag == 0){
      if(factor_flag == 0){
	printf("\n\n>>>>> LU method enabled >>>>>\n\n");
      }
      if(factor_flag == 1){
	printf("\n\n>>>>> Cholesky method enabled >>>>>\n\n");
      }
    }
    if(iter_flag == 1){
      if(factor_flag == 0){
	printf("\n\n>>>>> Bi-Cg method enabled >>>>>\n\n");
      }
      if(factor_flag == 1){
	printf("\n\n>>>>> Cg method enabled >>>>>\n\n");
      }
    }
  }
  if(sparse_flag == 1){
    if(iter_flag == 0){
      if(factor_flag == 0){
	printf("\n\n>>>>> LU method for sparse matrices enabled >>>>>\n\n");
      }
      if(factor_flag == 1){
	printf("\n\n>>>>> Cholesky method for sparse matrices enabled >>>>>\n\n");
      }
    }
    if(iter_flag == 1){
      if(factor_flag == 0){
	printf("\n\n>>>>> Bi-Cg method for sparse matrices enabled >>>>>\n\n");
      }
      if(factor_flag == 1){
	printf("\n\n>>>>> Cg method for sparse matrices enabled >>>>>\n\n");
      }
    }
  }
  
  /* ######################## STAGE 2 ######################## */
  //Create the MNA array and b vector
  if(sparse_flag == 0){
    mna(n, m2);
    if(tran_flag == 1){
      mna_tran();
    }
  }
  else if(sparse_flag == 1){
    mna_sparse(n, m2, n_z);
    if(tran_flag == 1){
      mna_sparse_tran(n, m2, C_n_z);
    }
  }
  /* ######################## PRINT ######################## */

  //print the list of the elements
  //printList();

  //print the hashtable with the nodes
  //printHashTable();

  printf("\n");

  printf("The number of nodes of DC analysis is: %d\n", n);
  printf("The number of elements in the second group is: %d\n", m2);
  //printf("The number of non zero elements of DC analysis is: %d\n", n_z);

  printf("\n");

  /*
  if(sparse_flag == 0){
    print_A_array();
    print_B_array();
    if(tran_flag == 1){
      	print_C_array();
      }
  }
 
  
  // ######################## STAGE 3 & STAGE 4 ######################## */
  //LU or Cholesky for non iterative methods  && Bi-CG or CG for iterative methods
  
  if((cholesky_flag == 1) && (factor_flag == 1) && (iter_flag == 0)){
    
    printf("The program must terminate!!!! You tried to run Cholesky algorithm in a non positive definite matrix...\n ");
    
    //free list
    free_list();

    //free hashtable
    freeHashTable();

    if(sparse_flag == 0){
      free(A_array);
      free(B_array);
      if(tran_flag == 1){
      	free(C_array);
      }
    }
    else if(sparse_flag == 1){
      free(B_array_cs);
      free(B_array_cs_temp);
    }
  
    freeList();

    freeList2();
  
    freeList3();
    
    return EXIT_FAILURE;
  
  }
  
  if(sparse_flag == 0){
    if(iter_flag == 0){	//flag for iterative methods
      //LU or Cholesky
      printf("\n\nRunning non iterative methods...\n\n ");

      factorization(factor_flag, n , m2);
      
      if(factor_flag == 0){
	lu_solve();
      }
      else if(factor_flag == 1){
	chol_solve();
      }
    
    
      /* ######################## PRINT ######################## */
      /*
      if(factor_flag == 0){
	printLU();
      }else{
	printChol();
      }
      
      */
      //printX();
    
    }
    else if(iter_flag == 1){	//flag for iterative methods

      printf("\n\nRunning iterative methods...\n\n ");
      
      if(factor_flag == 0){
	cg_method(sparse_flag, factor_flag, 0, itol, n, m2);
      }
      else if(factor_flag == 1){
	cg_method(sparse_flag, factor_flag, 0, itol, n, m2);
      }
      
      /* ######################## PRINT ######################## */
      
      //printIter_x(factor_flag);
      
    }
  }
  else if(sparse_flag == 1){
    if(iter_flag == 0){
      
      printf("\n\nRunning non iterative methods...\n\n ");
      
      clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);
      
      factorization_cs(factor_flag);
      
       if(factor_flag == 0){
	lu_solve_cs(0);
      }
      else if(factor_flag == 1){
	chol_solve_cs(0);
      }
      
       /* This is the end of the main computation. Take the end time,  *
	 * calculate the duration of the computation and report it. 	*/
      clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);
    
      printf ("\n\nTotal time in CPU for non iterative methods = %10g seconds\n\n",
		    (double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
		    (double) (tv2.tv_sec - tv1.tv_sec));
  
      /* ######################## PRINT ######################## */
      
     // printX_cs();
      
    }
    else if(iter_flag == 1){	//flag for iterative methods
      
      printf("\n\nRunning iterative methods...\n\n ");
      
      clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);
      
      if(factor_flag == 0){
	cg_method(sparse_flag, factor_flag, 0, itol, n, m2);
      }
      else if(factor_flag == 1){
	cg_method(sparse_flag, factor_flag, 0, itol, n, m2);
      }
      
      /* This is the end of the main computation. Take the end time,  *
	 * calculate the duration of the computation and report it. 	*/
      clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);
    
      printf ("\n\nTotal time in CPU for iterative methods = %10g seconds\n\n",
		    (double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
		    (double) (tv2.tv_sec - tv1.tv_sec));
      
      /* ######################## PRINT ######################## */
      
      printf("\n*****Sparse Matrices*****\n\n");
      
      //printIter_x(factor_flag);
    }
  }
  

     
  /* ################# dc sweep ################# */ 
  
  if(dc_sweep_flag == 1){ 
    m = 0;
    dc_tmp = root;
    plot_tmp = root2;
    //temp_B_array = tmp2_array;

    
    while(dc_tmp!=NULL){
      tmp = 0;
      tmp1 = 0;
      
      if(dc_tmp->element_id == 'v'){// dc sweep for Volt Source
	
	//find the k position of the Volt Source
	search_list(dc_tmp->input_variable);
	
	if(sparse_flag == 0){
	  if(iter_flag == 0){
	    tmp = B_array[n+sweep_k];
	    
	    //dc sweep for the given element
	    dc_sweep(dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  
	  }
	  else if(iter_flag == 1){
	    tmp = B_array[n+sweep_k];	//keep the n+sweep_k element of b vector
	  
	    iter_dc_sweep(dc_tmp->element_id, factor_flag, n, m2, dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  }
	}
	else if(sparse_flag == 1){
	  if(iter_flag == 0){
	    tmp = B_array_cs[n+sweep_k];
	    
	    //dc sweep for the given element
	    dc_sweep_cs(dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  
	  }
	  else if(iter_flag == 1){
	    
	    tmp = B_array_cs[n+sweep_k];	//keep the n+sweep_k element of b vector
	    
	    iter_dc_sweep_cs(dc_tmp->element_id, factor_flag, n, m2, dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  
	  }
	}
      }
      else if(dc_tmp->element_id == 'i'){// dc sweep for Power Source
	
	//find the positive and negative node of the Power Source
	search_list_nodes(dc_tmp->input_variable);
	
	if(sparse_flag == 0){
	  if(iter_flag == 0){
	    //dc sweep for the given element
	    dc_sweep2(dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value, factor_flag);
	  }
	  else if(iter_flag == 1){
	    //keep the original b vector elements in pos and neg position
	    if(pos != -1){
	      tmp = B_array[pos];
	    }
	    if(neg != -1){
	      tmp1 = B_array[neg];
	    }
	    
	    iter_dc_sweep(dc_tmp->element_id, factor_flag, n, m2, dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  }
	}
	else if(sparse_flag == 1){
	  if(iter_flag == 0){
	    //dc sweep for the given element
	    
	    dc_sweep2_cs(dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value, factor_flag);
	  }
	  else if(iter_flag == 1){
	    //keep the original b vector elements in pos and neg position
	    if(pos != -1){
	      tmp = B_array_cs[pos];
	    }
	    if(neg != -1){
	      tmp1 = B_array_cs[neg];
	    }
	    
	    iter_dc_sweep_cs(dc_tmp->element_id, factor_flag, n, m2, dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value);
	  }
	}
      }
      
      
      //print the DC sweep array
      //dc_sweep_printf(dc_tmp->input_variable, dc_tmp->element_id);
      
      
      if(plot_flag == 1){
	plot_option_print(dc_tmp->input_variable, dc_tmp->element_id, dc_tmp->start_value, dc_tmp->end_value, dc_tmp->inc_value, m);
      }
      
      m++;	//counter for the number of .DC - .PLOT
      
      //return the b vector in original state for another .DC - .PLOT
      if(dc_tmp->element_id == 'v'){
	if(sparse_flag == 0){
	  B_array[n + sweep_k] = tmp;
	}
	else if(sparse_flag == 1){
	  B_array_cs[n + sweep_k] = tmp;
	}
      }
      else if(dc_tmp->element_id == 'i'){
	if(sparse_flag == 0){
	  if(pos != -1){
	    B_array[pos] = tmp;
	  }
	  if(neg != -1){
	    B_array[neg] = tmp1;
	  }
	}
	if(sparse_flag == 1){
	  if(pos != -1){
	    B_array_cs[pos] = tmp;
	  }
	  if(neg != -1){
	    B_array_cs[neg] = tmp1;
	  }
	}
      }
      
      free_dc_sweep_array();	//free dc sweep array every time for different dc sweep
      
      dc_tmp = dc_tmp->next;	//take the  next element for dc sweep
    }
    
    
  }
  
  /* ######################## Stage 6 ######################## */
 
  if(tran_flag == 1){
    //if(sparse_flag == 1){
     // free(B_array_cs);
     // free(B_array_cs_temp);
     // free(x_vector);
     // cs_sfree(S);
     // cs_nfree(N);
   // }
    transient(n, m2, time_step, finish_time, method_flag, factor_flag, iter_flag, sparse_flag, itol);
  
    //tran_printf();
    tran_plot_print(time_step, finish_time);
  }
  
  /* ######################## Stage 7 ######################## */

  
  
  if(ac_flag == 1){
    
    printf("\n\n***** AC ANALYSIS *****\n\n");

    if(lin_flag == 1){
      double step;
      magnitudes = (double*)calloc(points*(n+m2), sizeof(double));
      phases = (double*)calloc(points*(n+m2), sizeof(double));
    
      if(sparse_flag == 0){
	if(iter_flag == 0){
	  if(factor_flag == 0){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    step = (end_freq )/ points;
	    
	    mna_ac(n, m2, start_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_lu();
	    
	    tmp_freq = start_freq + step;
	    
	    for(i=0; i < points; i++){
	      lu_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	     
	      //if(i==1){printX_ac();}
	      mna_ac(n, m2, tmp_freq);
	      
	      tmp_freq += step;
	    }
	  
	    ac_plot_print(start_freq, end_freq, points);
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = start_freq;
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_chol();
	    
	    
	    for(i=0; i < points; i++){
	      chol_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq += start_freq;
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printX_ac();
	}
	if(iter_flag == 1){
	  
	  if(factor_flag == 0){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	  
	    step = (end_freq )/ points;
	    
	    mna_ac(n, m2, start_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    tmp_freq = start_freq + step;
	    
	    for(i=0; i < points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	      tmp_freq += step;
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq, points);
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = start_freq;
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    print_A_ac_array();
	    print_B_ac_array();
	    
	    
	    for(i=0; i < points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    //abs arg polar
	    
	      tmp_freq += start_freq;
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      //if(i==1){printIter_x_ac(factor_flag);}
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printIter_x_ac(factor_flag);
	}
      }
      //for sparse matrices
      if(sparse_flag == 1){
	if(iter_flag == 0){
	  if(factor_flag == 0){
	    init_A_sparse(ac_n_z);	//initialise A and b complex arrays
	    init_B();
	    
	    step = (end_freq )/ points;
	    
	    mna_ac_sparse(n, m2, start_freq);	//create A and b 
	    
	    //print_A_B_array_cs_ac();
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    tmp_freq = start_freq + step;
	    
	    for(i=0; i < points; i++){
	      factorization_cs_ac(factor_flag);
	      
	      lu_solve_cs_ac(i);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(B_ac_array_temp[ii]);
		value_i = cimag(B_ac_array_temp[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    
	    //abs arg polar
	    
	      
	      
	      //compute the mna again for a different frequency
	      cs_ci_spfree(A_ac_sparse);
	      free(B_ac_array);
	      init_A_sparse(ac_n_z);
	      init_B();
	      
	      //if(i==1){printX_ac();}
	      mna_ac_sparse(n, m2, tmp_freq);
	      
	      tmp_freq += step;
	    }
	  
	    ac_plot_print(start_freq, end_freq, points);
	    
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = start_freq;
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_chol();
	    
	    
	    for(i=0; i < points; i++){
	      chol_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq += start_freq;
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printX_ac_cs();
	}
	if(iter_flag == 1){
	  
	  if(factor_flag == 0){
	    init_A_sparse(ac_n_z);	//initialise A and b complex arrays
	    init_B();
	  
	    step = (end_freq )/ points;
	    
	    mna_ac_sparse(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_B_array_cs_ac();
	    //print_A_ac_array();
	  // print_B_ac_array();
	    
	    tmp_freq = start_freq + step;
	    
	    for(i=0; i < points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      
	      
	      //compute the mna again for a different frequency
	      cs_ci_spfree(A_ac_sparse);
	      free(B_ac_array);
	      init_A_sparse(ac_n_z);
	      init_B();
	      
	      mna_ac_sparse(n, m2, tmp_freq);
	      
	      tmp_freq += step;
	    }
	  
	    ac_plot_print(start_freq, end_freq, points);
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = start_freq;
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    print_A_ac_array();
	    print_B_ac_array();
	    
	    
	    for(i=0; i < points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    //abs arg polar
	    
	      tmp_freq += start_freq;
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      //if(i==1){printIter_x_ac(factor_flag);}
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printIter_x_ac(factor_flag);
	}
      }
    }
    
    else if(log_flag == 1){
      
      magnitudes = (double*)calloc((points+1)*(n+m2), sizeof(double));
      phases = (double*)calloc((points+1)*(n+m2), sizeof(double));
    
      double vima;
     
      vima = (log10(end_freq) - log10(start_freq)) / points;
     
      
      if(sparse_flag == 0){
	if(iter_flag == 0){
	  if(factor_flag == 0){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(10, log10(start_freq));
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_lu();
	    
	    
	    for(i=0; i <= points; i++){
	      lu_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq = pow(10, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      //if(i==1){printX_ac();}
	      mna_ac(n, m2, tmp_freq);
	      
	   
	    }
	  
	    ac_plot_print_log(start_freq, end_freq, points);
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(start_freq, log10(start_freq) + vima);
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_chol();
	    
	    
	    for(i=0; i < points; i++){
	      chol_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq = pow(start_freq, log10(start_freq) + (vima * (i+2)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printX_ac();
	}
	if(iter_flag == 1){
	  
	  if(factor_flag == 0){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	  
	    tmp_freq = pow(10, log10(start_freq));
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    
	    
	    for(i=0; i <= points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq = pow(10, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	   
	   ac_plot_print_log(start_freq, end_freq, points);
	    
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(start_freq, log10(start_freq) + vima);
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    print_A_ac_array();
	    print_B_ac_array();
	    
	    
	    for(i=0; i < points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    //abs arg polar
	    
	       tmp_freq = pow(start_freq, log10(start_freq) + (vima * (i+2)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      //if(i==1){printIter_x_ac(factor_flag);}
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printIter_x_ac(factor_flag);
	}
      }
      //for sparse matrices
      if(sparse_flag == 1){
	if(iter_flag == 0){
	  if(factor_flag == 0){
	    init_A_sparse(ac_n_z);	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(10, log10(start_freq));
	    
	    mna_ac_sparse(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_B_array_cs_ac();
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    for(i=0; i <= points; i++){
	      factorization_cs_ac(factor_flag);
	      
	      lu_solve_cs_ac(i);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(B_ac_array_temp[ii]);
		value_i = cimag(B_ac_array_temp[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    
	    //abs arg polar
	    
	      tmp_freq = pow(10, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      cs_ci_spfree(A_ac_sparse);
	      free(B_ac_array);
	      init_A_sparse(ac_n_z);
	      init_B();
	      
	      //if(i==1){printX_ac();}
	      mna_ac_sparse(n, m2, tmp_freq);
	      
	      //printf("reeee\n");
	    }
	  
	    ac_plot_print_log(start_freq, end_freq, points);
	    
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(start_freq, log10(start_freq) + vima);
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    //print_A_ac_array();
	    //print_B_ac_array();
	    
	    initialize_chol();
	    
	    
	    for(i=1; i <= points; i++){
	      chol_fact();
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = GSL_REAL(gsl_vector_complex_get(x_ac, ii));
		value_i = GSL_IMAG(gsl_vector_complex_get(x_ac, ii));
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq = pow(start_freq, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printX_ac_cs();
	}
	if(iter_flag == 1){
	  
	  if(factor_flag == 0){
	    init_A_sparse(ac_n_z);	//initialise A and b complex arrays
	    init_B();
	  
	    tmp_freq = pow(10, log10(start_freq));
	    
	    mna_ac_sparse(n, m2, tmp_freq);	//create A and b 
	    
	   // print_A_B_array_cs_ac();
	    //print_A_ac_array();
	  // print_B_ac_array();
	    
	    
	    
	    for(i=0; i <= points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	      
	    //abs arg polar
	    
	      tmp_freq = pow(10, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      cs_ci_spfree(A_ac_sparse);
	      free(B_ac_array);
	      init_A_sparse(ac_n_z);
	      init_B();
	      
	      mna_ac_sparse(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print_log(start_freq, end_freq, points);
	  }
	  /*
	  if(factor_flag == 1){
	    init_A();	//initialise A and b complex arrays
	    init_B();
	    
	    tmp_freq = pow(start_freq, log10(start_freq) + vima);
	    
	    mna_ac(n, m2, tmp_freq);	//create A and b 
	    
	    print_A_ac_array();
	    print_B_ac_array();
	    
	    
	    for(i=1; i <= points; i++){
	      cg_method_ac(sparse_flag, factor_flag, i, itol, n, m2);
	      
	      //take the magnitudes and the phases from every step for the ac analysis plot
	      for(ii=0; ii<(n+m2); ii++){
		value_r = creal(iter_x_ac[ii]);
		value_i = cimag(iter_x_ac[ii]);
		magnitudes[j] = sqrt((value_r * value_r) + (value_i * value_i));
		phases[j] = atan(value_i/value_r) * 180 / PI;
		j++;
	      }
	    //abs arg polar
	    
	      tmp_freq = pow(start_freq, log10(start_freq) + (vima * (i+1)));
	      
	      //compute the mna again for a different frequency
	      free(A_ac_array);
	      free(B_ac_array);
	      init_A();
	      init_B();
	      
	      //if(i==1){printIter_x_ac(factor_flag);}
	      mna_ac(n, m2, tmp_freq);
	      
	    
	    }
	  
	    ac_plot_print(start_freq, end_freq);
	  }*/
	  
	  //printIter_x_ac(factor_flag);
	}
      }
    }
   // 
  
  }

  /* ######################## FREE MEMORY ######################## */

  //for(i=0; i < (points+1)*(n+m2); ++i){
    //printf("%g\n", magnitudes[i]);
  //}
  

  //free list
  free_list();

  //free hashtable
  freeHashTable();

  if(sparse_flag == 0){
    free(A_array);
    free(B_array);
    if(tran_flag == 1){
      	free(C_array);
	freeTran();
      }
  }
  else if(sparse_flag == 1){
    free(B_array_cs);
    free(B_array_cs_temp);
  }
    
  if(ac_flag == 1){
    if(sparse_flag = 0){
      free(A_ac_array);
    }
    else if(sparse_flag == 1){
      cs_ci_spfree(A_ac_sparse);
    }
    free(B_ac_array);
    free(magnitudes);
    free(phases);
  }
  //free(input_variable);
  
  if(sparse_flag == 0){
    if(iter_flag == 0){
      freeLU();
    }
    else if(iter_flag == 1){
      freeIter();
    }
  }
  else if(sparse_flag == 1){
    cs_spfree(Asparse);		//free the A sparse 
    cs_spfree(Csparse);		//free the C tilda sparse
    free(x_vector);
    cs_sfree(S);
    cs_nfree(N);
  }
  
  
  freeList();

  freeList2();
  
  freeList3();
  
  freeList4();

  if(ac_flag == 1){
    free(ac_method);
    //
    if(sparse_flag == 0){
      if(iter_flag == 0){
	if(factor_flag == 0){
	  freeLU_ac();
	}
	if(factor_flag == 1){
	  freeChol_ac();
	}
      }
      if(iter_flag == 1){
	
	freeIter_ac();
      }
    }
     if(sparse_flag == 1){
      if(iter_flag == 0){
	if(factor_flag == 0){
	  free(x_ac_cs);
	  free(B_ac_array_temp);
	}
	if(factor_flag == 1){
	  freeChol_ac();
	}
      }
      if(iter_flag == 1){
	
	freeIter_ac();
      }
    }
  }
  
  //free_dc_sweep_array();
  

  printf("\n*******************END OF PROGRAM************************\n");

  return EXIT_SUCCESS;
    
}
