double f(double *x, int n);
double f1(double x);

/*
int main(int argc, char** argv)
{
	double r = 10000000000.00;
	
	int n = 4;
	int m = 0;
	int p = 3;
	double* x = (double*) malloc( sizeof(double) * n );
    x[0]    = -0.01;
    x[1]    = +0.99;
	x[2]    = +1.99;
    x[3]    = -0.99;
	
	RnFunction *h = (RnFunction*) malloc(sizeof(RnFunction*) * m);
	RnFunction *g = (RnFunction*) malloc(sizeof(RnFunction*) * p);
	g[0] = g1;
	g[1] = g2;
	g[2] = g3;
	
	penalty_method(f, x, n, h, m, g, p, r);
	
	printf("%f %f %f %f\n", x[0], x[1], x[2], x[3]);
	
	free(x);
	free(h);
	free(g);
	
	return 0;
}

double f(double *x, int n)
{
	return x[0]*x[0] + x[1]*x[1] + 2*x[2]*x[2] + x[3]*x[3] - 5*x[0] - 5*x[1] - 21*x[2] + 7*x[3];
}

double g1(double *x, int n)
{
	return 8 - x[0]*x[0] - x[1]*x[1] - x[2]*x[2] -  x[3]*x[3] - x[0] + x[1] - x[2] + x[3];
}

double g2(double *x, int n)
{
	return 10 - x[0]*x[0] - 2*x[1]*x[1] - x[2]*x[2] - 2*x[3]*x[3] + x[0] - x[3];
}

double g3(double *x, int n)
{
	return 5 - 2*x[0]*x[0] - x[1]*x[1] - x[2]*x[2] - 2*x[0] + x[1] + x[3];
}
*/

int count = 0;
void gradient_test();
void straf_test();
void svenn_test();

double P(double *x, int n)
{
	return f(x,n);// + (1.0/r) * h(x,n) * h(x,n);
}

void test1(RnFunction *fs, double* x, int n)
{
	int i=0;
	for ( i = 0; i < n; i++)
	{
		printf("%f\n", fs[i](x, n));
	}
}

void lagrange(RnFunction f, double *x, int n, 
			  RnFunction* h, int k, double* u, 
			  RnFunction* g, int p, double* v, 
			  double* w);

void straf_test()
{
/*
    double epsilon	= 0.001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
	
	int n = 2;
	double* x = (double*)malloc( sizeof(double) * n );
	
    x[0]    = +5.0;
    x[1]    = +5.0;
	while ( r * (h(x,n)) > epsilon ) {
		conjugate_gradient_method(P, x, n, line_eps, gold_eps, grad_eps, epsilon);
		r = r * 0.1;
	}
	
	free(x);
*/
}

double f_rosenbrock(double *x, int n)
{
	double x1 = x[0];
	double x2 = x[1];
	count++;
    return ((1-x1)*(1-x1)) + 100*(x2-x1*x1)*(x2-x1*x1);
}

double f_1(double x)
{
	return (10.0 - x)*(10.0 - x);
}

double f_2(double x)
{
	return 2*x*x+16/x;
}

void svenn_test()
{
//	double dx = 0.5;
//	double x0 = 3.0;
//	double a = 6.0;
//	double b = 15.0;

//	search_interval_svenn(f_1, x0, dx, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	halph_interval_method(f_1, 0.001, &a, &b);
//	printf("a=%8.2f\nb=%8.2f\n", a, b);
//	search_method_pauella(f_2, 1.0, 1.0, 0.03, &a, &b);

	newton_raphson(f_2, 1.3333333, 0.0001);
}

void gradient_test()
{
    double epsilon	= 0.001;		//dovrun sona catma meyari
	double grad_eps	= 0.005;		//gradient
	double line_eps	= 0.1;			//parcani bolme
	double gold_eps	= 0.0001;		//qizil qayda ucun
    
	int n = 2;
    double* x  = (double*) malloc( sizeof(double) * n );
    x[0]    = -1.2;
    x[1]    = +1.0;
	conjugate_gradient_method(f_rosenbrock, x, n, line_eps, gold_eps, grad_eps, epsilon);
	free(x);
	
	printf("%d\n", count);
}

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}