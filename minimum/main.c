#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "methods.h"
#include "print.h"
#include "optimal.h"

extern void calculate();
extern void __calculate();

double func1(double *x, int n)
{
	return 2*x[0]*x[0] + x[1]*x[1] + 2*x[0]*x[1];
}

int main(int argc, char** argv)
{
	//sample_gradient1();
	__calculate();
	//calculate();
	/*
	double epsilon	= 0.1;		//dovrun sona catma meyari
	double grad_eps	= 0.00001;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.000001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    
	x[0]    = -6.0;
    x[1]    = -1.0;
	sample_project(func1, x, n, line_eps, gold_eps, grad_eps, epsilon, printer2);

	free(x);
	*/
	
	return 0;
}
