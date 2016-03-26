#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtable.h"

//flag for the first time, for the hashtable Initialization
int flag = 1;

//the node name into an integer
int id_num = 0;

int n=0;

/*Initialize hash table*/
void init_table(){
  int i;

  for(i = 0; i < SIZE; i++){
    hash_table[i] = NULL;
  }

}

unsigned long hash_matching(char *str){
  unsigned long hash = 0;
  int c;

  while (c = *str++){
    hash = c + (hash << 6) + (hash << 16) - hash;
  }
  
  return hash % SIZE;
}

void insert_hash_element(size_t keyID, char *str){
  
  struct hash_elements *curr_add;
  
  curr_add = (struct hash_elements *)malloc(sizeof(struct hash_elements));
  curr_add->key = (char *)malloc((strlen(str) + 1)*sizeof(char));
  
  //the first node of the list
  if(hash_table[keyID] == NULL){
    curr_add->next = NULL;
    hash_table[keyID] = curr_add; 
  }
  //put every node in the head of the list
  else{
    curr_add->next = hash_table[keyID];
    hash_table[keyID] = curr_add;
  }
  
  strcpy(curr_add->key, str);
  curr_add->id = id_num;
  
  //every node has its ordered int
  id_num++;


}


int search_element(size_t keyID, char *str){
  struct hash_elements *curr;

  //the hashtable is free
  if(hash_table[keyID] == NULL){
    return 1;
  }
  
  //the hashtable isn't free. Search if the node already exists
  for(curr = hash_table[keyID];  curr != NULL ; curr = curr->next){
    if(strcmp(curr->key,str) == 0){
      return 0;		//return 0 if the node exists
    }
  }
  
  //return 1 if nothing found
  return 1;

}

void insert(char *str, char c){
 //printf("-- %d\n", n);
  
  size_t keyID;
  
  if(flag){
    init_table();
    flag = 0;
  }
  
  if(strcmp(str, "0") == 0){
    return;
  }
  
  //hashing the node
  keyID = hash_matching(str);
  
  
  if(search_element(keyID, str)){
    //increase the counter for R,L,C,V,I elements nodes without ground node
    if((c != 'q') && (c != 'Q') && (c != 'm') && (c != 'M') && (c != 'd') && (c != 'D')){
      n++;
    }
    insert_hash_element(keyID, str);
  }

}


void freeHashTable(){
  int i;
  struct hash_elements *tmp;
  
  for(i=0; i<SIZE; i++){
    if(hash_table[i]!= NULL){
      while(hash_table[i] != NULL){
	tmp = hash_table[i];
	free(tmp->key);
	hash_table[i] = tmp->next;
	free(tmp);
      }
    }
  }
  
}



void printHashTable(){
  int i;
  struct hash_elements *tmp;
  
  for(i=0; i<SIZE; i++){
    if(hash_table[i]!= NULL){
      tmp = hash_table[i];
      printf("The #%d of the hashtable contains the nodes: ", i);
      while(tmp != NULL){
	printf("%s with id-> %d, ", tmp->key, tmp->id);
	tmp = tmp->next;
      }
      printf("\n");
    }
  }
}

int search_element_id(char *str){
  struct hash_elements *tmp;
  size_t keyID;
  
  keyID = hash_matching(str);
  
  tmp = hash_table[keyID];
  
  while(tmp!=NULL){
    if(strcmp(str, tmp->key) == 0){
      return tmp->id;	//return the id of the node
    }
    tmp = tmp->next;
  }
}