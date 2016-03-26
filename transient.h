#ifndef TRANSIENT_H
#define TRANSIENT_H


double *tran_array;

void transient(int n, int m2, double time_step, double finish_time, int m_flag, int f_flag, int i_flag, int s_flag, double itol);

void backward_euler(double time_step, double finish_time, int flag);

void backward_euler_iter(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol);

void backward_euler_iter_cs(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol);

void backward_euler_cs(double time_step, double finish_time, int flag);

void compute_E(double t);

void trapezoidal(double time_step, double finish_time, int flag);

void trapezoidal_cs(double time_step, double finish_time, int flag);

void trapezoidal_iter(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol);

void trapezoidal_iter_cs(int n, int m2, double time_step, double finish_time, int flag, int cs_flag, double itol);

double linear_value(double x1, double x2, double y1, double y2, double x);

void dgemv_tran(int flag, double a, double *array, double *xx, double *yy);

void dgemv_tran_cs(int flag, double a, cs *array, double *xx, double *yy);

void vector_add_tran(double a, double *xx, double *yy, double *zz);

void freeTran();

void tran_printf();

#endif