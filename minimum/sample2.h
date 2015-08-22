#ifndef _SAMPLE2_H_
#define _SAMPLE2_H_

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
	double **u1;
	
} Process2;

#endif //_SAMPLE2_H_
