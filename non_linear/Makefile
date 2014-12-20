CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g -Wall
CLDFLAGS	=    -Wl,-subsystem,console -mthreads
INCLUDES	=	-I..

all: main

main: main.o
	$(CC) $(CLDFLAGS) -o main.exe main.o

main.o: main.c
	$(CC) -c $(INCLUDES) -o main.o main.c
	
clean:
	del *.o *.exe

