CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g
INCLUDES	=	-I.
LDFLAGS		= 	-L. 

all: test_lib1 test_lib2 main

main: main.o
	$(CC) -o main.exe main.o $(LDFLAGS) -ltest_lib1

main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDES) -c main.c -o main.o
	
test_lib1: test_lib1.o
	$(CC) -Wl,-s -shared -mthreads -Wl,--out-implib,libtest_lib1.a -o test_lib1.dll test_lib1.o 

test_lib1.o: test_lib1.c
	$(CC) $(CFLAGS) $(INCLUDES) -c test_lib1.c -o test_lib1.o
	
test_lib2: test_lib2.o
	$(CC) -Wl,-s -shared -mthreads -Wl,--out-implib,libtest_lib2.a -o test_lib2.dll test_lib2.o 

test_lib2.o: test_lib2.c
	$(CC) $(CFLAGS) $(INCLUDES) -c test_lib2.c -o test_lib2.o
	
clean:
	del *.o *.a *.dll *.exe

