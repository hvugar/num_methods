CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g
INCLUDES	=	-I..

all: main

main: minimum.o gradient.o methods.o main.o
	$(CC) -o main.exe main.o minimum.o gradient.o methods.o $(LFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c -o main.o
	
minimum.o: minimum.c
	$(CC) $(CFLAGS) $(INCLUDES) -c minimum.c -o minimum.o
	
gradient.o: gradient.c
	$(CC) $(CFLAGS) $(INCLUDES) -c gradient.c -o gradient.o
	
methods.o: methods.c
	$(CC) $(CFLAGS) $(INCLUDES) -c methods.c -o methods.o
	
clean:
	del *.o *.exe

