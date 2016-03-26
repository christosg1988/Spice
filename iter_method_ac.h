#ifndef ITER_METHOD_AC_H
#define ITER_METHOD_AC_H

#include <complex.h>
#include "mna.h"
#include "csparse.h"
#include "mna_sparse.h"
#include "ac_analysis.h"
//#include "cs.h"

extern complex double *iter_x_ac;
complex double *M_ac;

void cg_method_ac(int cs_flag, int flag, int reset, double itol, int n, int m2);

void dgemv_ac(int flag, double a, complex double *array, complex double *xx, complex double *yy);

void dgemv_cs_ac(int flag, double a, cs_ci *array, complex double *xx, complex double *yy);

void vector_add_ac(int flag, complex double a, complex double *xx, complex double *yy, complex double *zz);

double norm_2_ac(complex double *xx);

complex double divide(complex double a, complex double b);


complex double ddot_ac(complex double *xx, complex double *yy);

void preconditioner_solve_ac(int flag, complex double *xx, complex double *yy, complex double *zz);

void init_cg_ac(int cs_flag, int flag);

void printIter_x_ac(int flag);

void freeIter_ac();

#endif