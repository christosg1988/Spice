#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "options.h"
#include "hashtable.h"
#include "dc_sweep.h"
#include "factorization.h"
#include "transient.h"
#include "ac_analysis.h"
#include "math.h"


dc *root = NULL;
dc *temp = NULL;

plot *root2 = NULL;
plot *temp2 = NULL;

t_plot *root3 = NULL;
t_plot *temp3 = NULL;

a_plot *root4 = NULL;
a_plot *temp4 = NULL;

double itol_t = 0.001;
int  factor_f = 0;
int sparse_f = 0;
int  iter_f = 0;
int method_f = 0;

dc* insert_node(dc *curr){
  
  curr = (dc*)malloc(sizeof(dc));
    
  //check for memory error
  if(curr == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  curr->next = NULL;
  
  //check if list is empty
  if(root == NULL){
    root = curr;
  }
  else{
    temp = root;
    
    while(temp->next != NULL){
      temp = temp->next;
    }
    
    temp->next = curr;
  }
  
  return(curr);
}

plot* insert_node2(plot *curr2){
  
  curr2 = (plot*)malloc(sizeof(plot));
    
  //check for memory error
  if(curr2 == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  curr2->next = NULL;
  
  //check if list is empty
  if(root2 == NULL){
    root2 = curr2;
  }
  else{
    temp2 = root2;
    
    while(temp2->next != NULL){
      temp2 = temp2->next;
    }
    
    temp2->next = curr2;
  }
  
  return(curr2);
}

void printDcList(){

  temp = root;
  
  while(temp!=NULL){
    printf("%s %g %g %g\n", temp->input_variable, temp-> start_value, temp->end_value, temp->inc_value);
  
    temp = temp->next;
  }
}

void printPlotList(){

  temp2 = root2;
  
  while(temp2!=NULL){
    printf("%s %d\n", temp2->plot_node, temp2->node_id);
  
    temp2 = temp2->next;
  }
}

void freeList(){
  
  while(root!=NULL){
    free(root->input_variable);
    
    temp = root;
    root = root->next;
    free(temp);
  }
}

void freeList2(){
  
  while(root2!=NULL){
    free(root2->plot_node);
    
    temp2 = root2;
    root2 = root2->next;
    free(temp2);
  }
}


void freeList3(){
  
  while(root3!=NULL){
    free(root3->tran_node);
    
    temp3 = root3;
    root3 = root3->next;
    free(temp3);
  }
}

void freeList4(){
  
  while(root4!=NULL){
    free(root4->ac_node);
    
    temp4 = root4;
    root4 = root4->next;
    free(temp4);
  }
}

/*Mporei na doulepsei kai me parapano apo ena Plot. Ypothetoume oti exoume
   osa dc sweep kai plot*/
void plot_option(char *token, int id){
  plot *curr2 = NULL;
  char del[] = " \t\n\r";
  
  while(token != NULL){
    
    if(token[0] == 'V' || token[0] == 'v' || token[0] == 'I' || token[0] == 'i'){
      
      curr2 = insert_node2(curr2);
      
      curr2->plot_node = (char*)malloc((strlen(token)-3)*sizeof(char));
      
      //take the name of the element for the plot
      strncpy(curr2->plot_node, &token[2], strlen(token) - 3);
      curr2->node_id = id;
      
    }
    
    token = strtok(NULL, del);
    
  }
  
}

/*Mporei na doulepsei kai me parapano apo ena Plot. Ypothetoume oti exoume
   osa tran kai plot*/
void tran_plot(char *token){
  char del[] = " \t\n\r()'V''v'";
  t_plot *new;
  
  token = strtok(NULL, del);
  
  while(token != NULL){
     new = (t_plot*)malloc(sizeof(t_plot));
     new->tran_node = (char*)malloc((strlen(token)+1)*sizeof(char));
     
     strcpy(new->tran_node, token);
     
     if(root3 == NULL){
       new->next = NULL;
       root3 = new;
    }
    else{
      new->next = root3;
      root3 = new;
    }
    
    token = strtok(NULL, del);
  }
  
}

/*Mporei na doulepsei kai me parapano apo ena Plot. Ypothetoume oti exoume
   osa ac kai plot*/
void ac_plot(char *token){
  char del[] = " \t\n\r()'V''v'";
  a_plot *new;
  
  token = strtok(NULL, del);
  
  while(token != NULL){
     new = (a_plot*)malloc(sizeof(a_plot));
     new->ac_node = (char*)malloc((strlen(token)+1)*sizeof(char));
     
     strcpy(new->ac_node, token);
     
     if(root4 == NULL){
       new->next = NULL;
       root4 = new;
    }
    else{
      new->next = root4;
      root4 = new;
    }
    
    token = strtok(NULL, del);
  }
  
}


/*ektyponei se diaforetika arxeia ola ta tran plot*/
void tran_plot_print(double step, double stop){
  int plot_node_id;
  int j, k, i;
  FILE *fp;
  char *name;
  
  t_plot *runner = root3;
  
  k = stop / step ;
  
  while(runner!=NULL){
    j = 0;
    
    plot_node_id = search_element_id(runner->tran_node);
    
    name = (char*)malloc((strlen(runner->tran_node)+6)*sizeof(char));
    
    //make the name of the output file
    name[0] = 'v';
    name[1] = '(';
    strcpy(&name[2], runner->tran_node);
    name[strlen(runner->tran_node) + 2] = ')';
    strcpy(&name[strlen(runner->tran_node)+3], ".gp");
    
    fp = fopen(name, "w+");	//we open a new file for every possible node
    
    if (fp == NULL){
      printf("Error opening file!\n");
      exit(1);
    }
    
    for(i=0; i <= k; i++){
      fprintf(fp, "%e\t%e\n", (i * step), tran_array[j+plot_node_id]);
      
      j += size;
    }
    
    fclose(fp);
    
    runner = runner->next;
    
    free(name);
  }

}

/*ektyponei se diaforetika arxeia ola ta ac plot gia grammiki sarosi*/
void ac_plot_print(double step, double stop, int points){
  int plot_node_id;
  int j, i;
  double k;
  FILE *fp;
  char *name;
  
  
  a_plot *runner = root4;
  
  k = (stop) / points ;
  
  while(runner!=NULL){
    j = 0;
    
    plot_node_id = search_element_id(runner->ac_node);
    
    name = (char*)malloc((strlen(runner->ac_node)+9)*sizeof(char));
    
    //make the name of the output file
    name[0] = 'v';
    name[1] = '(';
    strcpy(&name[2], runner->ac_node);
    name[strlen(runner->ac_node) + 2] = ')';
    strcpy(&name[strlen(runner->ac_node)+3], "_ac");
    strcpy(&name[strlen(runner->ac_node)+6], ".gp");
    
    fp = fopen(name, "w+");	//we open a new file for every possible node
    
    if (fp == NULL){
      printf("Error opening file!\n");
      exit(1);
    }
    
    for(i=0; i < points; i++){
      fprintf(fp, "%e\t%e\t%e\n", step + (i * k), magnitudes[j+plot_node_id], phases[j+plot_node_id]);
      
      j += size;
    }
    fclose(fp);
    
    runner = runner->next;
    
    free(name);
  }

}

/*ektyponei se diaforetika arxeia ola ta ac plot gia logarithmiki sarosi*/
void ac_plot_print_log(double step, double stop, int n){
  int plot_node_id;
  int j, i;
  double k;
  FILE *fp;
  char *name;
  
  a_plot *runner = root4;
  
  k = (log10(stop) - log10(step)) / n;
 
  while(runner!=NULL){
    j = 0;
    
    plot_node_id = search_element_id(runner->ac_node);
    
    name = (char*)malloc((strlen(runner->ac_node)+9)*sizeof(char));
    
    //make the name of the output file
    name[0] = 'v';
    name[1] = '(';
    strcpy(&name[2], runner->ac_node);
    name[strlen(runner->ac_node) + 2] = ')';
    strcpy(&name[strlen(runner->ac_node)+3], "_ac");
    strcpy(&name[strlen(runner->ac_node)+6], ".gp");
    
    fp = fopen(name, "w+");	//we open a new file for every possible node
    
    if (fp == NULL){
      printf("Error opening file!\n");
      exit(1);
    }
    
    
    for(i=0; i <= n; i++){
      fprintf(fp, "%e\t%e\t%e\n", (pow(10, log10(step) + (k * (i)))) , 20*log10(magnitudes[j+plot_node_id]), phases[j+plot_node_id]);
      
      j += size;
    }
    fclose(fp);
    
    runner = runner->next;
    
    free(name);
  }

}

/*ektyponei se diaforetika arxeia ola ta dc sweep*/
void plot_option_print(char *str, char id, double start, double stop, double inc, int m){
  int plot_node_id;
  double i;
  int j;
  FILE *fp;
  char *name;
  
  
  temp2 = root2;
  
  
  while(temp2!=NULL){
    j = 0;
    if(temp2->node_id == m){
      plot_node_id = search_element_id(temp2->plot_node);
      
      name = (char*)malloc((strlen(str)+strlen(temp2->plot_node)+5)*sizeof(char));
      
      //make the name of the output file
      if(id == 'v'){	//for Volt source element
	name[0] = 'v';
      }
      else if(id == 'i'){	//for Power source element
	name[0] = 'i';
      }
      
      strcpy(&name[1], str);
      name[strlen(str)+1] = '_';
      strcpy(&name[strlen(str)+2], temp2->plot_node);
      strcpy(&name[(strlen(str)+strlen(temp2->plot_node)+2)], ".gp");
      
      fp = fopen(name, "w+");	//we open a new file for every possible node
      
      if (fp == NULL){
	printf("Error opening file!\n");
	exit(1);
      }
      
      for(i=start; i < (stop+0.00000001); i += inc){
	fprintf(fp, "%g\t%g\n", i, dc_sweep_array[j+plot_node_id]);
	
	j += size;
      }
      fclose(fp);
    }
    
    
    temp2 = temp2->next;
    
  }

}

/*krataei se lista ola ta dc sweep pou mporei na exoume*/
int dc_option(char *token){
  int i = 0;
  dc *curr = NULL;
  char *value;
  int flag = 0;
  char del[] = " \t\n\r";
  
  if(token != NULL){
    curr = insert_node(curr);
    flag = 1;
  }
  
  while(token != NULL){
    
    if(token != NULL){
      if(i==0){
	//id for dc sweep for volt source and power source
	if(token[0] == 'V' || token[0] == 'v'){
	  curr->element_id = 'v';
	}
	else if(token[0] == 'I' || 'i'){
	  curr->element_id = 'i';
	}
	
	curr->input_variable = (char*)malloc((strlen(token)-1)*sizeof(char));
	
	strcpy(curr->input_variable, &token[1]);
	
      }
      else if(i==1){
	curr->start_value = strtod(token, &value);
      }
      else if(i==2){
	curr->end_value = strtod(token, &value);
      }
      else{
	curr->inc_value = strtod(token, &value);
	break;
      }
      i++;
    }
    token = strtok(NULL, del);
   }
   
   return flag;
}

void all_options(char *token){
  char *value;
  char del[] = " \t\n\r";
  
  while(token!=NULL){
    token = strtok(NULL, del);
    
    if(token!=NULL){
      if((strcmp(token, "SPD") == 0) || (strcmp(token, "spd") == 0)){
	factor_f = 1;
      }
      else if((strcmp(token, "ITER") == 0) || (strcmp(token, "iter") == 0)){
	iter_f = 1;
      }
      else if((strncmp(token, "ITOL", 4) == 0) || (strncmp(token, "itol", 4) == 0)){
	  itol_t = strtod(&token[5], &value);
      }
      else if((strcmp(token, "SPARSE") == 0) || (strcmp(token, "sparse") == 0)){
	sparse_f = 1;
      }
       else if((strcmp(token, "METHOD=TR") == 0) || (strcmp(token, "METHOD=BE") == 0)|| (strcmp(token, "method=tr") == 0)|| (strcmp(token, "method=be") == 0)){
	 if((strcmp(&token[7], "BE") == 0) || (strcmp(&token[7], "be") == 0)){
	  method_f = 1;
	}
      }
    }
    
  }
}