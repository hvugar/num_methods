CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g
INCLUDES	=	-I..

all: main

main: minimum.o main.o
	$(CC) -o main.exe main.o minimum.o $(LFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c -o main.o
	
minimum.o: minimum.c
	$(CC) $(CFLAGS) $(INCLUDES) -c minimum.c -o minimum.o
	

clean:
	del *.o *.exe

