#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "parser.h"
//#include "hashtable.h"

element *tmp = NULL;
element *head = NULL;
//int pwl_count = -1;
pwl *p_head = NULL;
pwl *p_new = NULL;
pwl *p_tmp = NULL;
  
element* insert_element(element *curr){
  
  curr = (element*)malloc(sizeof(element));
  
  //check for memory error
  if(curr == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  
  //check if list is empty
  if(head == NULL){
    curr->next = NULL;
    head = curr;
  }
  else{
    curr->next = head;
    head = curr;
  }
  
  return(curr);
}

void printList(){
  
  
  tmp = head;
  
  while(tmp!=NULL){
    
    //print the node id
    //printf("The node id is: %d\n", tmp->id);
    //print the node with Volt_Source structures
    if(tmp->type_id== 'v'){
      Volt_Source* temp = (Volt_Source*)(tmp->type);
      
      printf("*****-VOLT SOURCE ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the value of the element
      printf("The value number is: %g\n\n", temp->value);
      
      if(temp->t_flag == 1){
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  printf("PWL SPEC:(%g %g) ", t_temp->t1, t_temp->i1);
	  if(t_temp->p_root != NULL){
	    p_tmp = p_head;
	    while(p_tmp != NULL){
	      printf("(%g %g) ", p_tmp->t1, p_tmp->i1);
	      p_tmp = p_tmp->next;
	    } 
	   }
	  
	  printf("\n");
	}
	
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  printf("EXP SPEC:(%g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->i2, t_temp->td1, t_temp->tc1, t_temp->td2, t_temp->tc2);
	  printf("\n");
	}
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	 
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  printf("SIN SPEC:(%g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->ia, t_temp->fr, t_temp->td, t_temp->df, t_temp->ph);
	  printf("\n");
	}
	
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  printf("PULSE SPEC:(%g, %g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->i2, t_temp->td, t_temp->tr, t_temp->tf, t_temp->pw, t_temp->per);
	  printf("\n");
	}
      }
      printf("\nAC\n");
      //print the magnitude of the element
      printf("The magnitude is: %g\n", temp->ac_ptr->mag);
      //print the phase of the element
      printf("The phase is: %g\n\n", temp->ac_ptr->phase);
      
      printf("\n###############################\n\n");
    }

    //print the node with Power_Source structures
    else if(tmp->type_id== 'i'){
      Power_Source* temp = (Power_Source*)(tmp->type);
      
      printf("*****-POWER SOURCE ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the value of the element
      printf("The value number is: %g\n\n", temp->value);
      
      if(temp->t_flag == 1){
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  printf("PWL SPEC:(%g %g) ", t_temp->t1, t_temp->i1);
	  if(t_temp->p_root != NULL){
	    p_tmp = p_head;
	    while(p_tmp != NULL){
	      printf("(%g %g) ", p_tmp->t1, p_tmp->i1);
	      p_tmp = p_tmp->next;
	    } 
	   }
	  
	  printf("\n");
	}
	
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  printf("EXP SPEC:(%g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->i2, t_temp->td1, t_temp->tc1, t_temp->td2, t_temp->tc2);
	  printf("\n");
	}
	
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  printf("SIN SPEC:(%g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->ia, t_temp->fr, t_temp->td, t_temp->df, t_temp->ph);
	  printf("\n");
	}
	
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  printf("PULSE SPEC:(%g, %g, %g, %g, %g, %g, %g) ", t_temp->i1, t_temp->i2, t_temp->td, t_temp->tr, t_temp->tf, t_temp->pw, t_temp->per);
	  printf("\n");
	}
      }
      
      printf("\nAC\n");
      //print the magnitude of the element
      printf("The magnitude is: %g\n", temp->ac_ptr->mag);
      //print the phase of the element
      printf("The phase is: %g\n\n", temp->ac_ptr->phase);
      
      printf("\n###############################\n\n");
    }
    //free the node with Resistance structures
    else if(tmp->type_id== 'r'){
      Resistance* temp = (Resistance*)(tmp->type);
      
      printf("*****-RESISTANCE ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the value of the element
      printf("The value number is: %g\n\n", temp->value);
      
      printf("###############################\n\n");
    }
    //print the node with Capacity structures
    else if(tmp->type_id== 'c'){
      Capacity* temp = (Capacity*)(tmp->type);
      
      printf("*****-CAPACITY ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the value of the element
      printf("The value number is: %g\n\n", temp->value);
      
      printf("###############################\n\n");
    }
    //print the node with Inductance structures
    else if(tmp->type_id== 'l'){
      Inductance* temp = (Inductance*)(tmp->type);
      
      printf("*****-INDUCTANCE ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the value of the element
      printf("The value number is: %g\n\n", temp->value);
      
      printf("###############################\n\n");
    }
    //print the node with Diode structures
    else if(tmp->type_id== 'd'){
      Diode* temp = (Diode*)(tmp->type);
      
      printf("*****-DIODE ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the positive node of the element
      printf("The positive node is: %s\n", temp->pos_node);
      //print the negative node of the element
      printf("The negative node is: %s\n", temp->neg_node);
      //print the area of the element
      printf("The area number is: %g\n\n", temp->area);
      
      printf("###############################\n\n");
    }
    //print the node with MOS_transistor structures
    else if(tmp->type_id== 'm'){
      MOS_transistor* temp = (MOS_transistor*)(tmp->type);
      
      printf("*****-MOS TRANSISTOR ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the drain node of the element
      printf("The drain node is: %s\n", temp->drain);
      //print the gate node of the element
      printf("The gate node is: %s\n", temp->gate);
      //print the source node of the element
      printf("The source node is: %s\n", temp->source);
      //print the body of the element
      printf("The body node is: %s\n", temp->body);
      //print the length of the element
      printf("The length number is: %g\n", temp->length);
      //print the width of the element
      printf("The width number is: %g\n\n", temp->width);
      
      printf("###############################\n\n");
    }
    //print the node with BJT_transistor structures
    else if(tmp->type_id== 'q'){
      BJT_transistor* temp = (BJT_transistor*)(tmp->type);
      
      printf("*****-BJT TRANSISTOR ELEMENT-*****\n");
      //print the name of the element
      printf("The name of the element is: %s\n", temp->name);
      //print the collector node of the element
      printf("The collector node is: %s\n", temp->collector);
      //print the base node of the element
      printf("The base node is: %s\n", temp->base);
      //print the emitter node of the element
      printf("The emitter node is: %s\n", temp->emitter);
      //print area of the element
      printf("The area number is: %g\n\n", temp->area);
      
      printf("###############################\n\n");
    }
    tmp=tmp->next;
   }
}


void free_list(){
  
    
  while(head != NULL){
    
    //free the node with Volt_Source structures
    if(head->type_id== 'v'){
      Volt_Source* temp = (Volt_Source*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      if(temp->t_flag == 1){
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  if(t_temp->p_root == NULL){
	    free(t_temp);
	  }else{
	    while(t_temp->p_root != NULL){
	      p_tmp = p_head;
	      p_head = p_head->next;
	      free(p_tmp);
	      t_temp->p_root = p_head;
	    }
	    free(t_temp);
	  }
	}
	
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  free(t_temp);
	}
	
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  free(t_temp);
	}
	
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  free(t_temp);
	}
	
	free(temp->tran_type);
      }
      free(temp->ac_ptr);
      free(temp);
    }
    //free the node with Power_Source structures
    else if(head->type_id== 'i'){
      Power_Source* temp = (Power_Source*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      if(temp->t_flag == 1){
	if((strcmp(temp->tran_type, "pwl") == 0) || (strcmp(temp->tran_type, "PWL") == 0)){
	   
	  tran_pwl* t_temp = (tran_pwl*)(temp->t_type);
	  if(t_temp->p_root == NULL){
	    free(t_temp);
	  }else{
	    while(t_temp->p_root != NULL){
	      p_tmp = p_head;
	      p_head = p_head->next;
	      free(p_tmp);
	      t_temp->p_root = p_head;
	    }
	    free(t_temp);
	  }
	  
	}
	
	else if((temp->tran_type[0] == 'e') || (temp->tran_type[0] == 'E')){
	  tran_exp* t_temp = (tran_exp*)(temp->t_type);
	  free(t_temp);
	}
	
	else if((temp->tran_type[0] == 's') || (temp->tran_type[0] == 'S')){
	  tran_sin* t_temp = (tran_sin*)(temp->t_type);
	  free(t_temp);
	}
	
	else{
	  tran_pulse* t_temp = (tran_pulse*)(temp->t_type);
	  free(t_temp);
	}
	
	free(temp->tran_type);
      }
      free(temp->ac_ptr);
      free(temp);
    }
    //free the node with Resistance structures
    else if(head->type_id== 'r'){
      Resistance* temp = (Resistance*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      free(temp);
    }
    //free the node with Capacity structures
    else if(head->type_id== 'c'){
      Capacity* temp = (Capacity*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      free(temp);
    }
    //free the node with Inductance structures
    else if(head->type_id== 'l'){
      Inductance* temp = (Inductance*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      free(temp);
    }
    //free the node with Diode structures
    else if(head->type_id== 'd'){
      Diode* temp = (Diode*)(head->type);
      free(temp->name);
      free(temp->pos_node);
      free(temp->neg_node);
      free(temp);
    }
    //free the node with MOS_transistor structures
    else if(head->type_id== 'm'){
      MOS_transistor* temp = (MOS_transistor*)(head->type);
      free(temp->name);
      free(temp->drain);
      free(temp->gate);
      free(temp->source);
      free(temp->body);
      free(temp);
    }
    //free the node with BJT_transistor structures
    else if(head->type_id== 'q'){
      BJT_transistor* temp = (BJT_transistor*)(head->type);
      free(temp->name);
      free(temp->collector);
      free(temp->base);
      free(temp->emitter);
      free(temp);
    }
    
    tmp = head;
    head = head->next;
    free(tmp);
  }
}


void volt_info( Volt_Source *vSource, char *token, int *ground_flag){
  
  char del[] = " \r\t\n,()";
  int i=0;
  char *value, *ac_tmp=NULL;
  
  
  //ac_tmp = (char*)malloc(strlen(token)*sizeof(char));

  //initialize ac analysis variables to zero
  //malloc for ac struct
  vSource->ac_ptr = (ac*)malloc(sizeof(ac));
	
  //magnitude of vSource element
  vSource->ac_ptr->mag = 0;
  
  //phase for vSource element
  vSource->ac_ptr->phase = 0;
	  
	  
  //give the name of the Volt Source
  vSource->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(vSource->name == NULL){
    printf("Cannot allocate memory\n");
    exit(EXIT_FAILURE);
  }
  //name of the element
  strncpy(vSource->name, &token[1], strlen(token));  
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	vSource->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(vSource->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	 
	strcpy(vSource->pos_node, token);
	
	
	//check if ground node exists
	if(strcmp(vSource->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else if(i==1){
	//negative node of element
	vSource->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(vSource->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(vSource->neg_node, token);
	
	
	//check if ground node exists
	if(strcmp(vSource->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else if(i==2){
	//value of element
	vSource->value = strtod(token, &value);
      }
      else if(i==3){
	if((token[0] != 'e') && (token[0] != 'E') && (token[0] != 's') && (token[0] != 'S') && (token[0] != 'p') && (token[0] != 'P') ){
	  vSource->t_flag = 0;
	  
	  if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0)){
	  
	    //magnitude of vSource element
	    token = strtok(NULL, del);
	    vSource->ac_ptr->mag = strtod(token, &value);
	    
	    //phase for vSource element
	    token = strtok(NULL, del);
	    vSource->ac_ptr->phase = strtod(token, &value);
	    
	    break;
	  }
	
	  i++;
	  continue;
	}
	
	vSource->tran_type = (char*)malloc((strlen(token)+1)*sizeof(char));
	strcpy(vSource->tran_type, token);
	
	vSource->t_flag = 1;
	
	//in case that we have exp in Volt element
	if((strcmp(token, "exp") == 0) || (strcmp(token, "EXP") == 0)){
	  tran_exp *exp_el;
	  exp_el = (tran_exp*)malloc(sizeof(tran_exp));
	  vSource->t_type = exp_el;
	  
	  //i1
	  token = strtok(NULL, del);
	  exp_el->i1 = strtod(token, &value);
	  
	  //i2
	  token = strtok(NULL, del);
	  exp_el->i2 = strtod(token, &value);
	  
	  //td1
	  token = strtok(NULL, del);
	  exp_el->td1 = strtod(token, &value);
	  
	  //tc1
	  token = strtok(NULL, del);
	  exp_el->tc1 = strtod(token, &value);
	  
	  //td2
	  token = strtok(NULL, del);
	  exp_el->td2 = strtod(token, &value);
	  
	  //tc2
	  token = strtok(NULL, del);
	  exp_el->tc2 = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have sin in Volt element
	if((strcmp(token, "sin") == 0) || (strcmp(token, "SIN") == 0)){
	  tran_sin *sin_el;
	  sin_el = (tran_sin*)malloc(sizeof(tran_sin));
	  vSource->t_type = sin_el;
	  
	  //i1
	  token = strtok(NULL, del);
	  sin_el->i1 = strtod(token, &value);
	  
	  //ia
	  token = strtok(NULL, del);
	  sin_el->ia = strtod(token, &value);
	  
	  //fr
	  token = strtok(NULL, del);
	  sin_el->fr = strtod(token, &value);
	  
	  //td
	  token = strtok(NULL, del);
	  sin_el->td = strtod(token, &value);
	  
	  //df
	  token = strtok(NULL, del);
	  sin_el->df = strtod(token, &value);
	  
	  //ph
	  token = strtok(NULL, del);
	  sin_el->ph = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have pulse in Volt element
	if((strcmp(token, "pulse") == 0) || (strcmp(token, "PULSE") == 0)){
	 
	  tran_pulse *pulse_el;
	  pulse_el = (tran_pulse*)malloc(sizeof(tran_pulse));
	  vSource->t_type = pulse_el;
	   
	  //i1
	  token = strtok(NULL, del);
	  pulse_el->i1 = strtod(token, &value);
	  
	  //i2
	  token = strtok(NULL, del);
	  pulse_el->i2 = strtod(token, &value);
	  
	  //td
	  token = strtok(NULL, del);
	  pulse_el->td = strtod(token, &value);
	  
	  //tr
	  token = strtok(NULL, del);
	  pulse_el->tr = strtod(token, &value);
	  
	  //tf
	  token = strtok(NULL, del);
	  pulse_el->tf = strtod(token, &value);
	  
	  //pw
	  token = strtok(NULL, del);
	  pulse_el->pw = strtod(token, &value);
	  
	  //per
	  token = strtok(NULL, del);
	  pulse_el->per = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have pwl in Volt element
	if((strcmp(token, "pwl") == 0) || (strcmp(token, "PWL") == 0)){
	  //pwl_count++;
	  tran_pwl * pwl_el;
	  pwl_el = (tran_pwl*)malloc(sizeof(tran_pwl));
	  vSource->t_type = pwl_el;
	  
	  //t1
	  token = strtok(NULL, del);
	  pwl_el->t1 = strtod(token, &value);
	  
	  //i1
	  token = strtok(NULL, del);
	  pwl_el->i1 = strtod(token, &value);
	  
	  pwl_el->p_root = NULL;
	  
	  token = strtok(NULL, del);
	  
	  
	  while(token != NULL){
	    
	    if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0 )){
	      //i++;
	      ac_tmp = token;
	      break;
	    }
	    
	    p_new = (pwl*)malloc(sizeof(pwl));
	    
	    //check for memory error
	    if(p_new == NULL){
	      printf("Cannot allocate memory\n");
	      exit(EXIT_FAILURE);
	    }
	    
	    //t1
	    //token = strtok(NULL, del);
	    p_new->t1 = strtod(token, &value);
	    
	    //i1
	    token = strtok(NULL, del);
	    p_new->i1 = strtod(token, &value);
	    
	    if(p_head == NULL){
	      p_new->next = NULL;
	      p_head = p_new;
	      pwl_el->p_root = p_head;
	    }
	    else{
	      p_tmp = p_head;
	      while(p_tmp->next !=NULL){
		p_tmp = p_tmp->next;
	      }
	      p_tmp->next = p_new;
	      p_new->next = NULL;
	    }
	    
	    token = strtok(NULL, del);
	  }
	  
	  i++;
	  continue;
	}
      }
      else if(i==4){

	//for the case that we have pwl.check pwl func
	if(ac_tmp != NULL){
	  if((strcmp(ac_tmp, "ac") == 0 ) || (strcmp(ac_tmp, "AC") == 0 )){
	    //magnitude of vSource element
		  

	    vSource->ac_ptr->mag = strtod(token, &value);
	
	    //phase for vSource element
	    token = strtok(NULL, del);
	    vSource->ac_ptr->phase = strtod(token, &value);
	    
	    break;
	  }
	}
	
	if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0)){
	
	  //magnitude of vSource element
	  token = strtok(NULL, del);
	  vSource->ac_ptr->mag = strtod(token, &value);
	  
	  //phase for vSource element
	  token = strtok(NULL, del);
	  vSource->ac_ptr->phase = strtod(token, &value);
	  
	  break;
	}
      }
      i++;
    }
   
  }
  //free(ac_tmp);
}


void power_info( Power_Source *pSource, char *token, int *ground_flag){
  
  char del[] = " \r\t\n,()";
  int i=0;
  char *value, *ac_tmp=NULL;
  
  
  //ac_tmp = (char*)malloc(strlen(token)*sizeof(char));
  
  //initialize the ac analysis variables to zero
  //malloc for ac struct
  pSource->ac_ptr = (ac*)malloc(sizeof(ac));
	
  //magnitude of vSource element
  pSource->ac_ptr->mag = 0;
  
  //phase for vSource element
  pSource->ac_ptr->phase = 0;
	     
  //give the name of the Power Source
  pSource->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(pSource->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
  
  //name of the element
  strncpy(pSource->name, &token[1], strlen(token));
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	pSource->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(pSource->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(pSource->pos_node, token);
	
	//check if ground node exists
	if(strcmp(pSource->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//negative node of element
	pSource->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(pSource->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(pSource->neg_node, token);
	
	//check if ground node exists
	if(strcmp(pSource->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==2){
	//value of element
	pSource->value = strtod(token, &value);
      }
      else if(i==3){
	if((token[0] != 'e') && (token[0] != 'E') && (token[0] != 's') && (token[0] != 'S') && (token[0] != 'p') && (token[0] != 'P') ){
	  pSource->t_flag = 0;
	  
	   if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0 )){
	    
	    //magnitude of vSource element
	    token = strtok(NULL, del);
	    pSource->ac_ptr->mag = strtod(token, &value);
	    
	    //phase for vSource element
	    token = strtok(NULL, del);
	    pSource->ac_ptr->phase = strtod(token, &value);
	    
	    break;
	  }
	  
	  i++;
	  continue;
	}
	
	pSource->tran_type = (char*)malloc((strlen(token)+1)*sizeof(char));
	strcpy(pSource->tran_type, token);
	
	pSource->t_flag = 1;
	
	//in case that we have exp in Power element
	if((strcmp(token, "exp") == 0) || (strcmp(token, "EXP") == 0)){
	  tran_exp *exp_el;
	  exp_el = (tran_exp*)malloc(sizeof(tran_exp));
	  pSource->t_type = exp_el;
	  
	  //i1
	  token = strtok(NULL, del);
	  exp_el->i1 = strtod(token, &value);
	  
	  //i2
	  token = strtok(NULL, del);
	  exp_el->i2 = strtod(token, &value);
	  
	  //td1
	  token = strtok(NULL, del);
	  exp_el->td1 = strtod(token, &value);
	  
	  //tc1
	  token = strtok(NULL, del);
	  exp_el->tc1 = strtod(token, &value);
	  
	  //td2
	  token = strtok(NULL, del);
	  exp_el->td2 = strtod(token, &value);
	  
	  //tc2
	  token = strtok(NULL, del);
	  exp_el->tc2 = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have sin in Power element
	if((strcmp(token, "sin") == 0) || (strcmp(token, "SIN") == 0)){
	  tran_sin *sin_el;
	  sin_el = (tran_sin*)malloc(sizeof(tran_sin));
	  pSource->t_type = sin_el;
	  
	  //i1
	  token = strtok(NULL, del);
	  sin_el->i1 = strtod(token, &value);
	  
	  //ia
	  token = strtok(NULL, del);
	  sin_el->ia = strtod(token, &value);
	  
	  //fr
	  token = strtok(NULL, del);
	  sin_el->fr = strtod(token, &value);
	  
	  //td
	  token = strtok(NULL, del);
	  sin_el->td = strtod(token, &value);
	  
	  //df
	  token = strtok(NULL, del);
	  sin_el->df = strtod(token, &value);
	  
	  //ph
	  token = strtok(NULL, del);
	  sin_el->ph = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have pulse in Power element
	if((strcmp(token, "pulse") == 0) || (strcmp(token, "PULSE") == 0)){
	 
	  tran_pulse *pulse_el;
	  pulse_el = (tran_pulse*)malloc(sizeof(tran_pulse));
	  pSource->t_type = pulse_el;
	   
	  //i1
	  token = strtok(NULL, del);
	  pulse_el->i1 = strtod(token, &value);
	  
	  //i2
	  token = strtok(NULL, del);
	  pulse_el->i2 = strtod(token, &value);
	  
	  //td
	  token = strtok(NULL, del);
	  pulse_el->td = strtod(token, &value);
	  
	  //tr
	  token = strtok(NULL, del);
	  pulse_el->tr = strtod(token, &value);
	  
	  //tf
	  token = strtok(NULL, del);
	  pulse_el->tf = strtod(token, &value);
	  
	  //pw
	  token = strtok(NULL, del);
	  pulse_el->pw = strtod(token, &value);
	  
	  //per
	  token = strtok(NULL, del);
	  pulse_el->per = strtod(token, &value);
	  
	  i++;
	  continue;
	}
	
	//in case that we have pwl in Power element
	if((strcmp(token, "pwl") == 0) || (strcmp(token, "PWL") == 0)){
	  //pwl_count++;
	  
	  tran_pwl * pwl_el;
	  pwl_el = (tran_pwl*)malloc(sizeof(tran_pwl));
	  pSource->t_type = pwl_el;
	  
	  //t1
	  token = strtok(NULL, del);
	  pwl_el->t1 = strtod(token, &value);
	  
	  //i1
	  token = strtok(NULL, del);
	  pwl_el->i1 = strtod(token, &value);
	  
	  pwl_el->p_root = NULL;
	  
	  token = strtok(NULL, del);
	  
	  
	  while(token != NULL){
	    //printf("ok\n");
	    if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0 )){
	      //i++;
	      ac_tmp = token;
	      break;
	    }
	    
	    p_new = (pwl*)malloc(sizeof(pwl));
	    
	    //check for memory error
	    if(p_new == NULL){
	      printf("Cannot allocate memory\n");
	      exit(EXIT_FAILURE);
	    }
	    
	    //t1
	    //token = strtok(NULL, del);
	    p_new->t1 = strtod(token, &value);
	    
	    //i1
	    token = strtok(NULL, del);
	    p_new->i1 = strtod(token, &value);
	    
	    if(p_head == NULL){
	      //printf("ok1\n");
	      p_new->next = NULL;
	      p_head = p_new;
	      pwl_el->p_root = p_head;
	      
	    }
	    else{
	     
	      p_tmp = p_head;
	      while(p_tmp->next != NULL){
		p_tmp = p_tmp->next;
	      }
	      p_tmp->next = p_new;
	      p_new->next = NULL;
	      
	    }
	    
	    token = strtok(NULL, del);
	  }
	  
	  i++;
	  continue;
	}
      }
      else if(i==4){
		
	//for the case that we have pwl.check pwl func
	if(ac_tmp != NULL){
	  if((strcmp(ac_tmp, "ac") == 0 ) || (strcmp(ac_tmp, "AC") == 0 )){
	    
	    //magnitude of vSource element
	    pSource->ac_ptr->mag = strtod(token, &value);
      
	    //phase for vSource element
	    token = strtok(NULL, del);
	    pSource->ac_ptr->phase = strtod(token, &value);
	    break;
	  }
	} 
	
	if((strcmp(token, "ac") == 0 ) || (strcmp(token, "AC") == 0 )){
	   
	  //magnitude of vSource element
	  token = strtok(NULL, del);
	  pSource->ac_ptr->mag = strtod(token, &value);
	  
	  //phase for vSource element
	  token = strtok(NULL, del);
	  pSource->ac_ptr->phase = strtod(token, &value);
	  
	  break;
	}	
      }
      i++;
    }
   
  }
  //free(ac_tmp);
}


void resistance_info( Resistance *resistance, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *value;
  
  //give the name of the Resistance
  resistance->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(resistance->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
	
  strncpy(resistance->name, &token[1], strlen(token)-1);
  
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	resistance->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(resistance->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(resistance->pos_node, token);
	
	//check if ground node exists
	if(strcmp(resistance->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//negative node of element
	resistance->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(resistance->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(resistance->neg_node, token);
	
	//check if ground node exists
	if(strcmp(resistance->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else{
	//value of element
	resistance->value = strtod(token, &value);
	break;
      }
      i++;
    }
   
  }
}


void capacity_info( Capacity *capacity, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *value;
  
  //give the name of the Capacity
  capacity->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(capacity->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
	
  strncpy(capacity->name, &token[1], strlen(token)-1);
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	capacity->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(capacity->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(capacity->pos_node, token);
	
	//check if ground node exists
	if(strcmp(capacity->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//negative node of element
	capacity->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(capacity->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(capacity->neg_node, token);
	
	//check if ground node exists
	if(strcmp(capacity->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else{
	//value of element
	capacity->value = strtod(token, &value);
	break;
      }
      i++;
    }
   
  }
}


void inductance_info( Inductance *inductance, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *value;
  
  //give the name of the Inductance
  inductance->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(inductance->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
	
  strncpy(inductance->name, &token[1], strlen(token));
  
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	inductance->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(inductance->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(inductance->pos_node, token);
	
	//check if ground node exists
	if(strcmp(inductance->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//negative node of element
	inductance->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(inductance->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(inductance->neg_node, token);
	
	//check if ground node exists
	if(strcmp(inductance->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else{
	//value of element
	inductance->value = strtod(token, &value);
	break;
      }
      i++;
    }
   
  }
}



void diode_info( Diode *diode, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *value;
  
  //give the name of the Diode
  diode->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(diode->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
	
  strncpy(diode->name, &token[1], strlen(token));
  
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//positive node of element
	diode->pos_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(diode->pos_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(diode->pos_node, token);
	
	//check if ground node exists
	if(strcmp(diode->pos_node, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//negative node of element
	diode->neg_node = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(diode->neg_node == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	if(newLine(token) == -1){
	  strcpy(diode->neg_node, token);
	  
	    //check if ground node exists
	  if(strcmp(diode->neg_node, "0") == 0){
	    *ground_flag = 1;
	  }
	  
	  diode->area = 1;
	  
	  break;
	}
	
	strcpy(diode->neg_node, token);
	
	//check if ground node exists
	if(strcmp(diode->neg_node, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else{
	//area of element
	diode->area = strtod(token, &value);
	break;
      }
      i++;
    }
  }
}



void mos_info( MOS_transistor *mos_tran, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *length;
  char *width;
  
  //give the name of the MOS transistor
  mos_tran->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(mos_tran->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
  //the name of the element   
  strncpy(mos_tran->name, &token[1], strlen(token));
  
  
  while(token != NULL){
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//drain node of element
	mos_tran->drain = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(mos_tran->drain == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(mos_tran->drain, token);
	
	//check if ground node exists
	if(strcmp(mos_tran->drain, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//gate node of element
	mos_tran->gate = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(mos_tran->gate == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(mos_tran->gate, token);
	
	//check if ground node exists
	if(strcmp(mos_tran->gate, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else if(i==2){
	//source node of element
	mos_tran->source = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(mos_tran->source == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(mos_tran->source, token);
	
	//check if ground node exists
	if(strcmp(mos_tran->source, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==3){
	//source node of element
	mos_tran->body = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(mos_tran->body == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(mos_tran->body, token);
	
	//check if ground node exists
	if(strcmp(mos_tran->body, "0") == 0){
	  *ground_flag = 1;
	}
      }
      else if(i==4){
	//length of element
	mos_tran->length = strtod(&token[2], &length);
      }
      else{
	//width of element
	mos_tran->width = strtod(&token[2], &width);
	break;
      }
      i++;
    }
   
  }
  
}



void bjt_info( BJT_transistor *bjt_tran, char *token, int *ground_flag){
  
  char del[] = " \t\n";
  int i=0;
  char *value;
  
  //give the name of the BJT transistor
  bjt_tran->name = (char*)malloc((strlen(token)+1)*sizeof(char));
  
  //check for memory error
  if(bjt_tran->name == NULL){
      printf("Cannot allocate memory\n");
      exit(EXIT_FAILURE);
  }
  //name of the element
  strncpy(bjt_tran->name, &token[1], strlen(token));
  
  
  while(token != NULL){
    //printf("|%s|", token);
    token = strtok(NULL, del);
    if(token != NULL){
      if(i==0){
	//collector node of element
	bjt_tran->collector = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(bjt_tran->collector == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(bjt_tran->collector, token);
	
	//check if ground node exists
	if(strcmp(bjt_tran->collector, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==1){
	//base node of element
	bjt_tran->base = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(bjt_tran->base == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	strcpy(bjt_tran->base, token);
	
	//check if ground node exists
	if(strcmp(bjt_tran->base, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(i==2){
	//emitter node of element
	bjt_tran->emitter = (char*)malloc((strlen(token)+1)*sizeof(char));
	
	//check for memory error
	if(bjt_tran->emitter == NULL){
	   printf("Cannot allocate memory\n");
	   exit(EXIT_FAILURE);
	}
	
	if(newLine(token) == -1){
	  strcpy(bjt_tran->emitter, token);
	  
	  //check if ground node exists
	  if(strcmp(bjt_tran->emitter, "0") == 0){
	    *ground_flag = 1;
	  }
	  
	  bjt_tran->area = 1;
	  
	  break;
	}
	
	strcpy(bjt_tran->emitter, token);
	
	//check if ground node exists
	if(strcmp(bjt_tran->emitter, "0") == 0){
	  *ground_flag = 1;
	}
	
      }
      else if(token != NULL){
	//area of element
	bjt_tran->area = strtod(token, &value);
	break;
      }
      i++;
    }
   
  }
}

int newLine(char *str){
  int i;
  
  for(i = 0; i < strlen(str); i++){
    if((str[i] == '\n') || (str[i] == '\r')){
      str[i] = '\0';
      return -1;
    }
  }
  
  return 0;
}