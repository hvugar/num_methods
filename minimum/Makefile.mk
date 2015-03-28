CC          = gcc
CXX         = g++
DEPS        = main.c
OBJ         = main.o
CFLAGS      = -g -Wall
CLDFLAGS    = -Wl,-subsystem,console -mthreads
INCLUDES    = -I..
OBJ_DIR     = .
DEL         = del /f

all: main

main: print.o sample_functions.o sample_gradient.o sample_penalty.o method_penalty.o method_grad.o method_conj.o methods.o minimum.o main.o
	$(CC) $(CLDFLAGS) -o main.exe main.o minimum.o method_grad.o method_conj.o methods.o print.o method_penalty.o sample_functions.o sample_gradient.o sample_penalty.o $(LFLAGS)

main.o: main.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o main.o main.c
	
minimum.o: minimum.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o minimum.o minimum.c 
	
methods.o: methods.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o methods.o methods.c

method_conj.o: method_conj.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_conj.o method_conj.c 

method_grad.o: method_grad.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_grad.o method_grad.c

method_penalty.o: method_penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_penalty.o method_penalty.c
	
sample_penalty.o: sample_penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_penalty.o sample_penalty.c
	
sample_gradient.o: sample_gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_gradient.o sample_gradient.c
	
sample_functions.o: sample_functions.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_functions.o sample_functions.c

print.o: print.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o print.o print.c
	
clean:
	$(DEL) *.o
	$(DEL) *.exe

