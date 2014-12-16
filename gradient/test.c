double f(double *x, int n);
double f1(double x);

double f(double *x, int n)
{
    return pow(x[0],3) + 2*pow(x[1],2) - 3*x[0] - 4*x[1];
}

double f1(double x)
{
	return (x-0.133)*(x-0.133)+0.2;
}