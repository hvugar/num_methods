#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*R1Function)(double);
typedef double (*RnFunction)(double);

double f(double, double);

//Initial condition
double fi(double x);
//Boundary conditions
double m1(double t);
double m2(double t);

int main(int argc, char** argv)
{
	double dx = 0.10;
	double dt = 0.10;
	double a  = 1.00;
	
	int M = 10;
	int N = 10;
	
	double b[9] = {0.0};
	double A[9][9] = {{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0}};
	
	int j = 0;
	int i = 0;
	for (i=1; i<=M-1; i++)
	{
		if ( i == 1 )
			b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m1(dt*(j+1));
		else if ( i == (M-1) )
			b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1)) + (dt/(dx*dx))*m2(dt*(j+1));
		else
			b[i-1] = fi(dx*i) + dt*f(dx*i,dt*(j+1));
	}
	
//	for (i=0; i<M-1; i++)
//	{
//		printf("%8.4f\n", b[i]);
//	}
	
	for (i=1; i<=M-1; i++)
	{
		if ( i == 1 )
		{
			A[i-1][0] = 1 + (2*dt)/(dx*dx);
			A[i-1][1] = -dt/(dx*dx);
		} 
		else if ( i == (M-1) )
		{
			A[i-1][i-2] = -dt/(dx*dx);
			A[i-1][i-1] = 1 + (2*dt)/(dx*dx);
		}
		else
		{
			A[i-1][i-2] = -dt/(dx*dx);
			A[i-1][i-1] = 1 + (2*dt)/(dx*dx);
			A[i-1][i-0] = -dt/(dx*dx);
		}
	}
	
	for (i=0; i<M-1; i++)
	{
		for (j=0; j<M-1; j++)
			printf("%8.4f ", A[i][j]);
		printf("| %8.4f\n", b[i]);
	}
	
	return 0;
}

double f(double x, double t)
{
	return -1.0;
}

double fi(double x)
{
	return x*x;
}

//Boundary conditions
double m1(double t)
{
	return t;
}

double m2(double t)
{
	return 1.0+t;
}
