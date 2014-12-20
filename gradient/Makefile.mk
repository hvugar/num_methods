CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g -Wall
LDFLAGS		=   -Wl,-subsystem,console -mthreads
INCLUDES	=	-I..

all: main

main: print.o minimum.o gradient.o methods_grad.o methods_conj.o methods.o main.o
	$(CC) $(CLDFLAGS) -o main.exe main.o minimum.o gradient.o methods_grad.o methods_conj.o methods.o print.o $(LFLAGS)

main.o: main.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o main.o main.c
	
minimum.o: minimum.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o minimum.o minimum.c 
	
gradient.o: gradient.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o gradient.o gradient.c
	
methods.o: methods.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o methods.o methods.c

methods_conj.o: methods_conj.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o methods_conj.o methods_conj.c 

methods_grad.o: methods_grad.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o methods_grad.o methods_grad.c
	
print.o: print.c
	$(CC) -c $(CFLAGS) $(INCLUDES) -o print.o print.c

	
clean:
	del *.o *.exe

