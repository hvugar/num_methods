#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void penalty_sample1();
extern void penalty_sample2();
extern void penalty_sample4();
extern void sample_gradient1();

void func1(int n, double *x, double *f)
{
	*f = x[0]*x[0] + 2*x[1];
}

void constr1(int n, double *x, double *f)
{
	f[0] = x[0]*x[0] + 2*x[1];
	f[1] = x[0]*x[0] + 2*x[1];
}

void func2(int n, double *x, double *f)
{
	int m = 2;
	func1(n, x, f);
	
	double *h = malloc(sizeof(double)*m);
	constr1(n, x, h);
	
	int i;
	for (i=0; i<m; i++)
	{
		*f += h[i]*h[i];
	}
}

int main(int argc, char** argv)
{
	//penalty_sample4();
	sample_gradient1();
	
//	double *x = malloc(sizeof(double)*2);
//	x[0] =1.0;
//	x[1] = 2.0;
	
//	double f1=0;
//	func2(2, x, &f1);
//	printf("%f\n", f1);

	return 0;
}