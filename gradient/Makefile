CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g
INCLUDES	=	-I..

all: main

main: minimum.o gradient.o methods_grad.o methods_conj.o methods.o main.o
	$(CC) -o main.exe main.o minimum.o gradient.o methods.o $(LFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c -o main.o
	
minimum.o: minimum.c
	$(CC) $(CFLAGS) $(INCLUDES) -c minimum.c -o minimum.o
	
gradient.o: gradient.c
	$(CC) $(CFLAGS) $(INCLUDES) -c gradient.c -o gradient.o
	
methods.o: methods.c
	$(CC) $(CFLAGS) $(INCLUDES) -c methods.c -o methods.o

methods_conj.o: methods_conj.c
	$(CC) $(CFLAGS) $(INCLUDES) -c methods_conj.c -o methods_conj.o

methods_grad.o: methods_grad.c
	$(CC) $(CFLAGS) $(INCLUDES) -c methods_grad.c -o methods_grad.o

	
clean:
	del *.o *.exe

