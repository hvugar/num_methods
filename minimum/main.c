#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

extern void calculate();
extern void __calculate();

double func1(double x)
{
	return (x+2.0)*(x+2.0) + 1;
}

int main(int argc, char** argv)
{
	//sample_gradient1();
	__calculate();
	//calculate();
	return 0;
}
