#ifndef ITER_METHOD_H
#define ITER_METHOD_H

#include "mna.h"
#include "csparse.h"
#include "mna_sparse.h"

extern double *iter_x;
double *M;

void cg_method(int cs_flag, int flag, int reset, double itol, int n, int m2);

void dgemv(int flag, double a, double *array, double *xx, double *yy);

void dgemv_cs(int flag, double a, cs *array, double *xx, double *yy);

void vector_add(double a, double *xx, double *yy, double *zz);

double norm_2(double *xx);

double ddot(double *xx, double *yy);

void preconditioner_solve(double *xx, double *yy, double *zz);

void init_cg(int cs_flag, int flag);

void printIter_x(int flag);

void freeIter();

#endif