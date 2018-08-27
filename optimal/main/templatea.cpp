#include "templatea.h"

void TemplateA::Main(int argc, char **argv)
{
    TemplateA tmp;
    DoubleVector y;
    tmp.calculate(0.001, 1000, y, 0.0, 0.0);

    FILE* file = fopen("d:/data.txt", "w");
    IPrinter::printVector(y, NULL, 1000, 0, 0, file);
    fclose(file);

    //DoubleVector y0(11);
    //for (unsigned int n=0; n<=10; n++) y0[n] = tmp.fx(n*0.1);
    //IPrinter::printVector(y0);
}

double TemplateA::fx(double x) const
{
    return x*x*x*x - 0.2*x*x*x - 0.8*x*x + 0.1;
}

double TemplateA::r(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 1.0; }

double TemplateA::p(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 0.0; }

double TemplateA::q(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const { return 0.0; }

double TemplateA::f(unsigned int i UNUSED_PARAM, double x UNUSED_PARAM) const
{
    if (i==800)
        return -5.0*1000.0;
    else
        return 0.0;

    //double sigma = 0.001;
    //return -5.0*(1.0/(sqrt(2.0*M_PI)*sigma)) * exp(-((x-0.8)*(x-0.8))/(2.0*sigma*sigma));
    //return r(i, x)*(12.0*x*x - 1.2*x - 1.6) +
    //       p(i,x)*(4.0*x*x*x - 0.6*x*x - 1.6*x) +
    //       q(i,x)*(x*x*x*x - 0.2*x*x*x - 0.8*x*x + 0.1);
}

void TemplateA::calculate(double h, unsigned int N, DoubleVector &y, double A, double B) const
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
