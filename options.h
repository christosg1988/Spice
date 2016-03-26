#ifndef OPTIONS_H
#define OPTIONS_H

struct list_opt{
  struct list_opt *next;
  char *input_variable;
  char element_id;
  double start_value;
  double end_value;
  double inc_value;
};

typedef struct list_opt dc;

struct list_opt2{
  struct list_opt2 *next;
  char *plot_node;
  int node_id;
};

typedef struct list_opt2 plot;

struct list_opt3{
  struct list_opt3 *next;
  char* tran_node;
};

typedef struct list_opt3 t_plot;

struct list_opt4{
  struct list_opt4 *next;
  char* ac_node;
};

typedef struct list_opt4 a_plot;

extern dc *root;

extern plot *root2;

extern t_plot *root3;

extern a_plot *root4;

extern int factor_f;
extern int iter_f;
extern int sparse_f;
extern double itol_t;
extern int method_f;

dc* insert_node(dc *curr);

double *magnitudes, *phases;

void freeList();

int dc_option(char *token);

void tran_plot(char *token);

void ac_plot(char *token);

void tran_plot_print(double step, double stop);

void ac_plot_print(double step, double stop, int points);

void freeList3();

void freeList4();

void all_options(char *token);

plot* insert_node2(plot *curr2);

void freeList2();

void ac_plot_print_log(double step, double stop, int n);

void plot_option(char *token, int id);

void plot_option_print(char *str, char id, double start, double stop, double inc, int m);

void printDcList();

void printPlotList();

#endif