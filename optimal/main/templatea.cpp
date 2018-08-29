#include "templatea.h"

unsigned int N = 100000;
double hx = 0.00001;
void TemplateA::Main(int argc, char **argv)
{
    TemplateA tmp;
    DoubleVector y;
    //tmp.calculateC(0.01, 100, y, 0.0, 0.0);
    tmp.calculateC(hx, N, y, tmp.fx(0.0), tmp.fd(0.0));
    IPrinter::printVector(y);

    FILE* file = fopen("d:/data.txt", "w");
    IPrinter::printVector(y, NULL, N, 0, 0, file);
    fclose(file);

    DoubleVector y0(11);
    for (unsigned int n=0; n<=10; n++) y0[n] = tmp.fx(n*0.1);
    IPrinter::printVector(y0);
}

double TemplateA::fx(double x) const
{
//    return x*x*x*x - 0.2*x*x*x - 0.8*x*x + 0.1;
    return exp(2.0*x-2.0) - 0.1;
}

double TemplateA::fd(double x) const
{
//    return x*x*x*x - 0.2*x*x*x - 0.8*x*x + 0.1;
    return 2.0*exp(2.0*x-2.0);
}

double TemplateA::r(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 1.0; }

double TemplateA::p(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 0.5; }

double TemplateA::q(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 0.3; }

double TemplateA::f(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const
{
    double fx = exp(2.0*x-2.0)*(4.0*r(i,x)+2.0*p(i,x)+q(i,x)) - 0.1*q(i,x);
    return fx;

    //if (i==60000) return fx += 5.0*N; else return fx += 0.0;

    //double sigma = h;
    //return fx += 5.0*(1.0/(sqrt(2.0*M_PI)*sigma)) * exp(-((x-0.6)*(x-0.6))/(2.0*sigma*sigma));

    //return r(i, x)*(12.0*x*x - 1.2*x - 1.6) +
    //       p(i,x)*(4.0*x*x*x - 0.6*x*x - 1.6*x) +
    //       q(i,x)*(x*x*x*x - 0.2*x*x*x - 0.8*x*x + 0.1);
}

void TemplateA::calculateB(double h, unsigned int N, DoubleVector &y, double A, double B) const
{
    double *a1 = (double*) malloc(sizeof(double)*(N-1));
    double *b1 = (double*) malloc(sizeof(double)*(N-1));
    double *c1 = (double*) malloc(sizeof(double)*(N-1));
    double *d1 = (double*) malloc(sizeof(double)*(N-1));
    double *rx = (double*) malloc(sizeof(double)*(N-1));

    for (unsigned int n=1; n<=N-1; n++)
    {
        double x = n*h;

        a1[n-1] = r(n,x) - 0.5*p(n,x)*h;
        b1[n-1] = -2.0*r(n,x) + q(n,x)*h*h;
        c1[n-1] = r(n,x) + 0.5*p(n,x)*h;
        d1[n-1] = f(n,x)*h*h;
    }

    d1[0]   -= a1[0]*A;
    d1[N-2] -= c1[N-2]*B;
    a1[0] = c1[N-2] = 0.0;

    tomasAlgorithm(a1, b1, c1, d1, rx, N-1);
    y.clear();
    y.resize(N+1);
    for (unsigned int i=0; i<=N-2; i++) y[i+1] = rx[i];
    y[0] = A; y[N] = B;

    free(rx);
    free(d1);
    free(c1);
    free(b1);
    free(a1);
}

void TemplateA::calculateC(double h, unsigned int N, DoubleVector &y, double a, double b) const
{
    y.resize(N+1);
    double x0 = 0.0;
    double x1 = h;

    y[0] = a;
    y[1] = a + h*b + ((h*h)/2.0) * ((p(0, x0)/r(0, x0))*b + (q(0, x0)/r(0, x0))*a + (f(0, x0)/r(0, x0)));
    y[1] += 5.0*(1.0/h)*(h*h*0.5);

    for (unsigned int i=2; i<=N; i++)
    {
        double x = i*h;
        double k = 1.0/(1.0+0.5*h*(p(i,x)/r(i,x)));
        y[i] = (2.0 - (q(i,x)/r(i,x))*h*h)*k*y[i-1] + (0.5*(p(i,x)/r(i,x))*h - 1.0)*k*y[i-2] + (f(i,x)/r(i,x))*h*h*k;
    }
}
