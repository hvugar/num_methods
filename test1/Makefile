CC			=	gcc
CXX			=	g++
DEPS		=	main.c
OBJ			=	main.o
CFLAGS		=	-g
INCLUDES	=	-I..

all: test_lib1 test_lib2 main

main: main.o
	$(CC) -o main.exe main.o $(LFLAGS)

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


#methods.o: methods.c
#	$(CC) $(CFLAGS) $(INCLUDES) -c methods.c -o methods.o

#methods_conj.o: methods_conj.c
#	$(CC) $(CFLAGS) $(INCLUDES) -c methods_conj.c -o methods_conj.o

#methods_grad.o: methods_grad.c
#	$(CC) $(CFLAGS) $(INCLUDES) -c methods_grad.c -o methods_grad.o

	
clean:
	del *.o *.exe

