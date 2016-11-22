CC          = gcc
CXX         = g++
CFLAGS      = -O2 -Wall
CXXFLAGS    = -O2 -Wall
LDFLAGS     = -L. 
INCLUDES    = -I.

all: main 

main: projection.o doublevector.o r1minimize.o parabolicequation.o function.o printer.o gradient.o gradient_cjt.o example31.o example32.o example33.o example34.o example35.o example36.o example37.o example38.o example335.o example336.o main.o
	$(CXX) -g -o main.exe projection.o doublevector.o r1minimize.o parabolicequation.o function.o printer.o gradient.o gradient_cjt.o example31.o example32.o example33.o example34.o example35.o example36.o example37.o example38.o example335.o example336.o main.o $(LDFLAGS)

main.o: main.cpp
	$(CXX) -c $(CXXFLAGS) -o main.o main.cpp

example31.o: example31.cpp
	$(CXX) -c $(CXXFLAGS) -o example31.o example31.cpp
	
example32.o: example32.cpp
	$(CXX) -c $(CXXFLAGS) -o example32.o example32.cpp
	
example33.o: example33.cpp
	$(CXX) -c $(CXXFLAGS) -o example33.o example33.cpp
	
example34.o: example34.cpp
	$(CXX) -c $(CXXFLAGS) -o example34.o example34.cpp
	
example35.o: example35.cpp
	$(CXX) -c $(CXXFLAGS) -o example35.o example35.cpp
	
example36.o: example36.cpp
	$(CXX) -c $(CXXFLAGS) -o example36.o example36.cpp
	
example37.o: example37.cpp
	$(CXX) -c $(CXXFLAGS) -o example37.o example37.cpp
	
example38.o: example38.cpp
	$(CXX) -c $(CXXFLAGS) -o example38.o example38.cpp
	
example335.o: example335.cpp
	$(CXX) -c $(CXXFLAGS) -o example335.o example335.cpp
	
example336.o: example336.cpp
	$(CXX) -c $(CXXFLAGS) -o example336.o example336.cpp
	
projection.o: projection.cpp
	$(CXX) -c $(CXXFLAGS) -o projection.o projection.cpp

doublevector.o: doublevector.cpp
	$(CXX) -c $(CXXFLAGS) -o doublevector.o doublevector.cpp
	
r1minimize.o: r1minimize.cpp
	$(CXX) -c $(CXXFLAGS) -o r1minimize.o r1minimize.cpp
	
function.o: function.cpp
	$(CXX) -c $(CXXFLAGS) -o function.o function.cpp
	
parabolicequation.o: parabolicequation.cpp
	$(CXX) -c $(CXXFLAGS) -o parabolicequation.o parabolicequation.cpp
		
printer.o: printer.cpp
	$(CXX) -c $(CXXFLAGS) -o printer.o printer.cpp
		
gradient.o: gradient.cpp
	$(CXX) -c $(CXXFLAGS) -o gradient.o gradient.cpp
	
gradient_cjt.o: gradient_cjt.cpp
	$(CXX) -c $(CXXFLAGS) -o gradient_cjt.o gradient_cjt.cpp	
	
clean:
	del *.o *.txt *.exe *.a
	
