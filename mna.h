#ifndef MNA_H
#define MNA_H

double* A_array;
double* B_array;
double* C_array;


void init_A_array();
void init_C_array();
void init_B_array();
void print_A_array();
void print_C_array();
void print_B_array();
void mna(int n, int m2);
void mna_tran();


#endif