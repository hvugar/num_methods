#ifndef _OPTIMAL4_H_
#define _OPTIMAL4_H_

#include <stdio.h>

typedef struct
{
	double dx;
	double dt;
	
	unsigned int n;

	double *x;
	double *t;
	
	double x0;
	double x1;
	
	double t0;
	double t1;
	
	double *f;
	
	double T1;
	
	double epsilon;
	
	double *p;
} Process4;

#endif // _OPTIMAL4_H_