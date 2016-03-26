#ifndef DC_SWEEP_H
#define DC_SWEEP_H

#include "factorization.h"
#include "mna.h"
#include "hashtable.h"
#include "parser.h"
#include "iter_method.h"
#include "options.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>
#include "mna_sparse.h"

//array for the outputs of x vector
double *dc_sweep_array;

//extern double sweep_size;
int sweep_k;

int pos;
int neg;


void dc_sweep(double start_value, double end_value, double inc_value);	//sweep for Volt source with LU method

void iter_dc_sweep(char t, int flag, int n, int m2, double start_value, double end_value, double inc_value);	//sweep with iterative methods using Cg or Bi-Cg methods

void iter_dc_sweep_cs(char t, int flag, int n, int m2, double start_value, double end_value, double inc_value);	//sweep with iterative methods using Cg or Bi-Cg methods from csparse

void dc_sweep2(double start_value, double end_value, double inc_value, int flag);	//sweep for Power source

void dc_sweep2_cs(double start_value, double end_value, double inc_value, int flag);	//sweep for Power source

void dc_sweep_cs(double start_value, double end_value, double inc_value);	//sweep for Volt source with Lu method from csparse

void dc_sweep_printf(char *str, char c);

void init_dc_sweep_array();

void free_dc_sweep_array();

void search_list(char *str);

void search_list_nodes(char *str);

#endif