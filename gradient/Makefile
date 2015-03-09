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

main: sample_functions.o sample_gradient.o sample_penalty.o penalty.o print.o minimum.o gradient.o methods_grad.o methods_conj.o methods.o main.o
	$(CC) $(CLDFLAGS) -o main.exe $(OBJ_DIR)/main.o $(OBJ_DIR)/minimum.o $(OBJ_DIR)/gradient.o $(OBJ_DIR)/methods_grad.o $(OBJ_DIR)/methods_conj.o $(OBJ_DIR)/methods.o $(OBJ_DIR)/print.o $(OBJ_DIR)/penalty.o $(OBJ_DIR)/sample_functions.o $(OBJ_DIR)/sample_gradient.o $(OBJ_DIR)/sample_penalty.o $(LFLAGS)

main.o: main.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/main.o main.c
	
minimum.o: minimum.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/minimum.o minimum.c 
	
gradient.o: gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/gradient.o gradient.c
	
methods.o: methods.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/methods.o methods.c

methods_conj.o: methods_conj.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/methods_conj.o methods_conj.c 

methods_grad.o: methods_grad.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/methods_grad.o methods_grad.c
	
print.o: print.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/print.o print.c
	
penalty.o: penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/penalty.o penalty.c
	
sample_penalty.o: sample_penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/sample_penalty.o sample_penalty.c
	
sample_gradient.o: sample_gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/sample_gradient.o sample_gradient.c
	
sample_functions.o: sample_functions.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o $(OBJ_DIR)/sample_functions.o sample_functions.c
	
clean:
	$(DEL) $(OBJ_DIR)\*.o
	$(DEL) *.exe

