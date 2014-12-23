CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g -Wall
CLDFLAGS	=   -Wl,-subsystem,console -mthreads
INCLUDES	=	-I..

all: main

main: diff_method.o tomas.o main.o
	$(CC) $(CLDFLAGS) -o main.exe main.o tomas.o diff_method.o

main.o: main.c
	$(CC) -c $(INCLUDES) -o main.o main.c
	
tomas.o: tomas.c
	$(CC) -c $(INCLUDES) -o tomas.o tomas.c	

diff_method.o: diff_method.c
	$(CC) -c $(INCLUDES) -o diff_method.o diff_method.c	
	
clean:
	del *.o *.exe

