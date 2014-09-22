CC		=	gcc
CXX		=	g++
DEPS	=	main.c
OBJ		=	main.o
CFLAGS	=	-g
INCLUDES=	-I..

all: main

main: main.o
	$(CC) -o main.exe main.o $(LFLAGS)

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c -o main.o

clean:
	del *.o main.exe main.s

asm:
	$(CC) $(INCLUDES) -S -masm=intel -m32 -fno-asynchronous-unwind-tables -o main.s main.c

