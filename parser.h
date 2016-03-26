#ifndef PARSER_H
#define PARSER_H


//struct for the list
struct list_el{
  struct list_el *next;
  void *type;  
  char type_id;
  //int id;
};

typedef struct list_el element;

extern element *head;


typedef struct{
  double mag;
  double phase;
}ac;


//struct that describes an indepedent Volt source
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  double value;
  int t_flag; //flag for the existence of transient type in the element(if flag==0 the element has its dc value for the transient analysis)
  char *tran_type;
  void *t_type;
  ac *ac_ptr;
}Volt_Source;



//struct that describes an indepedent Power source
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  double value;
  int t_flag; //flag for the existence of transient type in the element(if flag==0 the element has its dc value for the transient analysis)
  char *tran_type;
  void *t_type;
  ac *ac_ptr;
}Power_Source;


//Struct that describes resistance
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  double value;
}Resistance;


//Struct that describes capacity
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  double value;
}Capacity;


//Struct that describes inductance
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  double value;
}Inductance;


//Struct that describes diode
typedef struct{
  char *name;
  char *pos_node;
  char *neg_node;
  //char *model_name;
  double area;
}Diode;


//Struct that describes MOS transistor
typedef struct{
  char *name;
  char *drain;
  char *gate;
  char *source;
  char *body;
  //char *model_name;
  double length;
  double width;
}MOS_transistor;


//Struct that describes BJT transistor
typedef struct{
  char *name;
  char *collector;
  char *base;
  char *emitter;
  //char *model_name;
  double area;
  
}BJT_transistor;

typedef struct{
  double i1;
  double i2;
  double td1;
  double tc1;
  double td2;
  double tc2;
}tran_exp;


typedef struct{
  double i1;
  double ia;
  double fr;
  double td;
  double df;
  double ph;
}tran_sin;


typedef struct{
  double i1;
  double i2;
  double td;
  double tr;
  double tf;
  double pw;
  double per;
}tran_pulse;




//struct for the list for pwl
struct list_pwl{
  struct list_pwl *next;
  double t1, i1;
};

typedef struct list_pwl pwl;

//struct for pwl
typedef struct{
  double t1, i1;
  pwl *p_root;
}tran_pwl;

extern pwl *p_head;

element* insert_element(element *curr);
void printList();
void free_list();
void volt_info(Volt_Source *vSource, char *token, int *ground_flag);
void power_info(Power_Source *pSource, char *token, int *ground_flag);
void resistance_info(Resistance *resistance, char *token, int *ground_flag);
void capacity_info(Capacity *capacity, char *token, int *ground_flag);
void inductance_info(Inductance *inductance, char *token, int *ground_flag);
void diode_info(Diode *diode, char *token, int *ground_flag);
void mos_info(MOS_transistor *mos_tran, char *token, int *ground_flag);
void bjt_info(BJT_transistor *bjt_tran, char *token, int *ground_flag);
int newLine(char *str);


#endif