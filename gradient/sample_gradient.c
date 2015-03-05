#include "methods.h"
#include "print.h"

extern double f_rosenbrock(double *x, int n);

void sample_gradient1()
{
    double epsilon	= 0.0001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    
	x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon, printer1);
	
	puts("");
	
    //x[0]    = -1.2;
    //x[1]    = +1.0;
	//conjugate_gradient_method1(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);

	free(x);
}
