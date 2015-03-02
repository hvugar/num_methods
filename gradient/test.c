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

double f1(double x)
{
	return (x-0.133)*(x-0.133)+0.2;
}

double straight_line_search_metod1(R1Function fx, double x0, double dx, double *a, double *b, int *count)
{
	double y0 = 0.0;
	double y1 = 0.0;
	double y2 = 0.0;
	
	// if at next and last point of x0 function is greater
	// then decrease the dx to half
	y0 = fx( x0 );
	(*count)++;
	y1 = fx( x0 - dx );
	(*count)++;
	y2 = fx( x0 + dx );
	(*count)++;
	
	while (y1 > y0 && y2 > y0)
	{
		dx = dx / 2.0;
		y1 = fx( x0 - dx );
		y2 = fx( x0 + dx );
		(*count)++;
		(*count)++;
	}
	
	if ( y2 > y0 )
	{
		dx = -1 * dx;
		y2 = fx( x0 + dx );
		(*count)++;

	}

    double x1 = x0;
    double x2 = x1 + dx;
	
    y1 = fx( x1 );
	(*count)++;
    y2 = fx( x2 );
	(*count)++;
	
    int i = 1;
    while ( y2 <= y1 )
    {
		x1 = x2;
        y1 = y2;
		
		i++;
        x2 = x0 + i * dx;
        y2 = fx(x2);

		(*count)++;
		
		dx = dx*2;
    }
	
	if ( dx < 0 )
	{
		*a = x0 + (i) * dx;
		*b = x0 + (i-2) * dx;	
	}
	
	
	if ( dx > 0 )
	{
	    *a = x0 + (i-2) * dx;
		*b = x0 + (i+0) * dx;
	}
	
	return (*a+*b)/2;
}