#ifndef _OPTIMAL1_H_
#define _OPTIMAL1_H_

#include "method_grid.h"

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
	
	double **f;
	
} Process1;

#endif //_OPTIMAL1_H_