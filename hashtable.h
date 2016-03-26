#ifndef HASHTABLE_H
#define HASHTABLE_H

#define SIZE 670494
//#define SIZE 1670494

struct hash_elements
{
	unsigned int id;
	char *key;
	struct hash_elements *next;
};

struct hash_elements *hash_table[SIZE];


//counter for R,L,C,V,I element nodes
extern int n;

void init_table();
unsigned long hash_matching(char *str);
void insert_hash_element(size_t keyID, char *str);
//void delete_element(unsigned char *str, int  scope);
int serch_element(size_t keyID, char *str);
void insert(char *str, char c);
void freeHashTable();
void printHashTable();
int search_element_id(char *str);




#endif