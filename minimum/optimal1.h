#ifndef _OPTIMAL1_H_
#define _OPTIMAL1_H_

#include "method_grid.h"
#include "print.h"

typedef struct {
	double *x;
	double *t;
	
	double x0;
	double x1;
	
	double t0;
	double t1;
	
	double dx;
	double dt;
	
	double alpha;
	int n;
	int m;
	
	double **f;
	double **u;
	double **p;
	double **g;
	
} Process1;

#endif //_OPTIMAL1_H_