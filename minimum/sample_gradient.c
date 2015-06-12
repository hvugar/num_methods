#include "methods.h"
#include "print.h"

extern double f_rosenbrock(double *x, int n);

void sample_gradient1()
{
    double epsilon	= 0.000001;		//dovrun sona catma meyari
	double grad_eps	= 0.000001;		//gradient
	double step     = 0.1;			//parcani bolme
	double gold_eps	= 0.000001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    
	x[0]    = -1.0;
    x[1]    = +1.2;
	//steepest_descent_gradient_method(f_rosenbrock, x, n, step, gold_eps, grad_eps, epsilon, printer1);
	conjugate_gradient_method(f_rosenbrock, x, n, step, gold_eps, grad_eps, epsilon, printer2);

	free(x);
}
