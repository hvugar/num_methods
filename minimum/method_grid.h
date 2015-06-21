#ifndef _GRID_METHOD_H_
#define _GRID_METHOD_H_

#include "methods.h"

typedef struct 
{
	double **u;
	int n;
	int m;
} grid;

void implicit_difference_scheme(R2Function f, R1Function fi, R1Function m1, R1Function m2, double alpha, double dx, double dt, double x1, double x2, double t1, double t2, grid *g);

#endif //_GRID_METHOD_H_