#include "optimal1.h"

double y(double x)
{
	return x*x;
}

double u(double x, double t)
{
	return x*x*x*x + t*t*t - 1.0;
}

double JSum() 
{
	
}