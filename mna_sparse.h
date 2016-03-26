#ifndef MNA_SPARSE_H
#define MNA_SPARSE_H
#include "csparse.h"

cs *Asparse;
cs *Csparse;
cs *C_comp;
cs *C_tilda_comp;
css *S;
csn *N;

double* B_array_cs;
double* B_array_cs_temp;
double* x_vector;


void mna_sparse(int n, int m2, int n_z);

void mna_sparse_tran(int n, int m2, int n_z);

void factorization_cs(int flag);

void lu_solve_cs(int flag);

void chol_solve_cs(int flag);

void print_A_B_array_cs();

void printX_cs();

#endif