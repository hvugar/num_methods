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
	double a, b, c;
	a = b = 0.0;
	c = straight_line_search_metod(func1, 0, 0.1, &a, &b);
	
	printf("%f %f %f\n", a, c, b);
	__calculate();
	//calculate();
	return 0;
}
