# Makefile. If you change it, remember than in makefiles multiple spaces
# ARE NOT EQUIVALENT to tabs. The line after a rule starts with a tab!

#Add any executable you want to be created here.
EXECUTABLES = myspice

#This is the compiler to use
CC = gcc

#These are the flags passed to the compiler. Change accordingly
GCCFLAGS = -g -I /usr/local/include -I ./CXSparse/Include -L ./CXSparse/Lib/ 

# make all will create all executables
all: $(EXECUTABLES) 

# This is the rule to create any executable from the corresponding .c 
# file with the same name.
all:	build

build:	myspice

main.o:	main.c
	$(CC) $(GCCFLAGS) -o $@ -c $< 
	
parser_func.o:	parser_func.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

hashtable.o: hashtable.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

mna.o: mna.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

factorization.o: factorization.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

dc_sweep.o: dc_sweep.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

options.o: options.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

iter_method.o:	iter_method.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

csparse.o:	csparse.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

mna_sparse.o:	mna_sparse.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

transient.o:	transient.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

ac_analysis.o: ac_analysis.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

iter_method_ac.o: iter_method_ac.c
	$(CC) $(GCCFLAGS) -o $@ -c $<

myspice:	main.o	parser_func.o hashtable.o mna.o factorization.o dc_sweep.o options.o iter_method.o mna_sparse.o csparse.o transient.o iter_method_ac.o ac_analysis.o
	$(CC) $(GCCFLAGS) -o $@ $+ -lgsl -lgslcblas -lm -lrt -lcxsparse




#this is a rule to run the executable
run:
	-./$(EXECUTABLES) $(INPUT)

# make clean will remove all executables 
clean:
	 rm -f $(EXECUTABLES) *.o *.gp 
