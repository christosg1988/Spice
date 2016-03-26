#ifndef AC_ANALYSIS_H
#define AC_ANALYSIS_H
#include <complex.h>
#include "cs.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex_math.h>

complex double *A_ac_array;
complex double *B_ac_array;
complex double *B_ac_array_temp;

cs_ci *A_ac_sparse;
cs_ci *C_ac_comp;
cs_cis *S_ac;
cs_cin *N_ac;

//cs_cl* Asparse_ac;
gsl_vector_complex *x_ac;

complex double *x_ac_cs;

void init_A();

void init_B();

//void mna_ac_cs(int n, int m2, double freq, int n_z);

void mna_ac(int n, int m2, double freq);

void print_A_ac_array();

void print_B_ac_array();

void freeLU_ac();

void printX_ac();

void initialize_lu();

void factorization_cs_ac(int flag);

void lu_solve_cs_ac(int flag);

void lu_fact();

void chol_fact();

void printX_ac_cs();

void initialize_chol();

void freeChol_ac();

void print_A_B_array_cs_ac();

void init_A_sparse(int n_z);

void mna_ac_sparse(int n, int m2, double freq);
#endif