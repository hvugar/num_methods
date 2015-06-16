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

dll: print.o sample_functions.o sample_gradient.o sample_penalty.o runga_kutta.o optimal.o optimal2.o method_penalty.o method_grad.o method_conj.o methods.o minimum.o main.o
	$(CC) -shared -o minimum.dll -Wl,--out-implib=libminimum.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import minimum.o method_grad.o method_conj.o methods.o print.o method_penalty.o runga_kutta.o optimal.o optimal2.o sample_functions.o sample_gradient.o sample_penalty.o

main: method_prj_grad.o print.o sample_functions.o sample_gradient.o sample_penalty.o runga_kutta.o optimal.o optimal2.o method_penalty.o method_grad.o method_conj.o methods.o minimum.o main.o
	$(CC) $(CLDFLAGS) -o main.exe main.o minimum.o method_grad.o method_conj.o methods.o print.o method_penalty.o runga_kutta.o optimal.o optimal2.o sample_functions.o sample_gradient.o sample_penalty.o method_prj_grad.o $(LFLAGS)

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

runga_kutta.o: runga_kutta.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o runga_kutta.o runga_kutta.c

optimal.o: optimal.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o optimal.o optimal.c
	
optimal2.o: optimal2.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o optimal2.o optimal2.c
	
sample_penalty.o: sample_penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_penalty.o sample_penalty.c
	
sample_gradient.o: sample_gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_gradient.o sample_gradient.c
	
sample_functions.o: sample_functions.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_functions.o sample_functions.c
	
sample_project.o: sample_project.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_project.o sample_project.c
	
method_prj_grad.o: method_prj_grad.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_prj_grad.o method_prj_grad.c

print.o: print.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o print.o print.c
	
clean:
	$(DEL) *.o
	$(DEL) *.exe
	$(DEL) *.dll
	$(DEL) *.dll.a
