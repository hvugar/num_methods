CC          = gcc
CXX         = g++
DEPS        = main.c
OBJ         = main.o
CFLAGS      = -g -Wall
CLDFLAGS    = -Wl,-subsystem,console -mthreads
INCLUDES    = -I..
OBJ_DIR     = .
DEL         = del /f
SMPL_OBJ    = sample_functions.o sample_gradient.o sample_penalty.o
OBJECTS     = method_prj_grad.o print.o $(SMPL_OBJ) runga_kutta.o optimal.o method_penalty.o method_grad.o method_conj.o method_grid.o methods.o minimum.o main.o

all: main

dll: $(OBJECTS)
	$(CC) -shared -o minimum.dll -Wl,--out-implib=libminimum.dll.a -Wl,--export-all-symbols -Wl,--enable-auto-import $(OBJECTS)

main: $(OBJECTS)
	$(CC) $(CLDFLAGS) -o main.exe $(OBJECTS) $(LFLAGS)

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
	
method_grid.o: method_grid.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_grid.o method_grid.c

optimal.o: optimal1.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o optimal.o optimal1.c
	
sample_penalty.o: sample_penalty.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_penalty.o sample_penalty.c
	
sample_gradient.o: sample_gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_gradient.o sample_gradient.c
	
sample_functions.o: sample_functions.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_functions.o sample_functions.c
	
sample_grid.o: sample_grid.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o sample_grid.o sample_grid.c
	
method_prj_grad.o: method_prj_grad.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o method_prj_grad.o method_prj_grad.c

print.o: print.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o print.o print.c
	
clean:
	$(DEL) *.o
	$(DEL) *.exe
	$(DEL) *.dll
	$(DEL) *.dll.a
