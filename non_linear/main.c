#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*R1Function)(double);

double f1(double x);
double g1(double x);
double h1(double x);

double newton_method(R1Function f, R1Function g, double a, double b, double x0, double epsilon)
{
	int k = 0;
	double x1 = x0 - f(x0)/g(x0);
	
	printf("%4d %8.4f %8.4f %8.4f %8.4f\n", k, x0, f(x0), g(x0), -f(x0)/g(x0));
	while ( fabs(x1 - x0) > epsilon )
	{
		x0 = x1;
		x1 = x0 - f(x0)/g(x0);
		k++;
		printf("%4d %8.4f %8.4f %8.4f %8.4f\n", k, x0, f(x0), g(x0), -f(x0)/g(x0));
	}
	return x1;
}

double half_div_method(R1Function f, double a, double b, double x0, double epsilon)
{
	int k = 0;
	double c = 0.0;
	double f_a = 0.0;
	double f_b = 0.0;
	double f_c = 0.0;
	
	printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", k, a, b, f(a), f(b), (a + b) / 2.0, f((a + b) / 2.0));
	
	while ( fabs(a - b) > (2.0 * epsilon) )
	{
		c = (a + b) / 2.0;
		double f_a = f(a);
		double f_b = f(b);
		double f_c = f(c);
		
		if (f_a * f_c > 0.0)
			a = c;
		if (f_b * f_c > 0.0)
			b = c;
			
		k++;
		printf("%4d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", k, a, b, f(a), f(b), (a + b) / 2.0, f((a + b) / 2.0));
	}
	
	return (a+b) / 2;
}

int main(int argc, char** argv)
{
	double a = 0.4;
	double b = 0.6;
	double epsilon = 0.001;
	double x0 = 0.6;
	double x = newton_method(f1, g1, a, b, x0, epsilon);
	puts("-----------------------");
	printf("%8.6f %8.6f\n", x, f1(x));
	return 0;
}

double f1(double x)
{
	return pow(M_E, 2*x) + 3*x - 4.0;
}

double g1(double x)
{
	return 2*pow(M_E, 2*x)+3;
}

double h1(double x)
{
	return 4*pow(M_E, 2*x);
}

double f(double x)
{
	return 0.5*x*sin(x) + 1.0;
}

double g(double x)
{
	return 0.5*(sin(x)+x*cos(x));
}

double h(double x)
{
	return 0.5*(2*cos(x)-x*sin(x));
}