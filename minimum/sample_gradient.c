#include "methods.h"
#include "print.h"
#include "gradient.h"

extern double f_rosenbrock(double *x, int n);

void sample_gradient1()
{
    double epsilon1  = 0.000001;   //dovrun sona catma meyari
    double epsilon2  = 0.000001;   //dovrun sona catma meyari
    double grad_eps  = 0.000001;   //gradient
    double line_step = 0.1;        //parcani bolme
    double gold_eps  = 0.000001;   //qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    
	x[0]    = -1.0;
    x[1]    = +1.2;

    SteepestDescentMethod(f_rosenbrock, Gradient, x, n, line_step, gold_eps, grad_eps, epsilon1, epsilon2);

	free(x);
}
