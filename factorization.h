#ifndef FACTORIZATION_H
#define FACTORIZATION_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_permutation.h>

extern gsl_vector *x;
extern gsl_permutation * p;
extern double *tmp_array;
extern double *tmp2_array;

extern int size;

void init_LU();

void init_Chol();

void factorization(int flag, int n, int m2);

void lu_Decomp();

void lu_solve();

void chol_Decomp();

void chol_solve();

void printLU();

void printChol();

void printX();

void freeLU();

#endif