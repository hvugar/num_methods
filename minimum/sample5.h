#ifndef _SAMPLE5_H_
#define _SAMPLE5_H_

#include <stdio.h>

typedef struct
{
	double dx;
	double dt;
	
	unsigned int N;

	double *t;
	double *x;
	double *psi;
	
	double t0;
	double t1;

	double x0;
	double x1;
	
	double psi0;
	double psi1;
	
	
	double T[2];
	double p[2];
	double grad[2];

	double epsilon;
} Process5;

#endif // _SAMPLE5_H_
