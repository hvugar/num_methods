#include "example7.h"

void Example7::Main(int argc, char *argv[])
{
    Example7 e;
    e.calculateK4();
}

Example7::Example7()
{
    h = 0.01;
    N = 100;
    K = 4;
    w = 14;
    p = 10;
}

void Example7::calculateK4()
{
    DoubleVector rx(N+1);
    for (unsigned int k=0; k<=N; k++) rx[k] = f(k);
    IPrinter::printVector(w,p,rx);

    DoubleMatrix betta(N+1, N+1);
    DoubleVector eta(N+1);

    betta[0][0] = 1.0;
    betta[0][N/2] = 3.5;
    betta[0][N] = -2.5;
    eta[0] = 0.0; for (unsigned int k=0; k<=N; k++) eta[0] += betta[0][k]*rx[k];
    printf("%.10f\n", eta[0]);

    for (unsigned int k=1; k<=N; k++)
    {

    }
}

double Example7::a(unsigned int k) const
{
    //double t = k*h;
    return 2.0;
}

double Example7::b(unsigned int k) const
{
    double t = k*h;
    return -2.0*sin(10.*t) - 10.0*cos(10.0*t);
}

double Example7::f(unsigned int k) const
{
    double t = k*h;
    return sin(10.0*t);
}
