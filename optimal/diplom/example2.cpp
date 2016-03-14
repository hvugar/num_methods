#include "example2.h"

// u(x,t)=sin(x)+e^t

void BorderParabolic1D2::main()
{
    DoubleMatrix m;
    BorderParabolic1D2 bp;
    bp.calculateU(m, bp.hx, bp.ht, bp.N, bp.M, 2.0);
    IPrinter::printMatrix(m);
    FILE *file = fopen("example2.txt", "w");
    IPrinter::printVector(m[m.size()-1], NULL, bp.N, 0, 0, file);
    fclose(file);
}

BorderParabolic1D2::BorderParabolic1D2()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);
    e1 = 0.20;
    e2 = 0.55;
    e3 = 0.84;
}

double BorderParabolic1D2::fi(unsigned int i) const
{
    double x = i*hx;
    return sin(x)+1.0;
}

double BorderParabolic1D2::m1(unsigned int j) const
{
    double t = j*ht;
    return exp(t);
}

double BorderParabolic1D2::m2(unsigned int j) const
{
    double t = j*ht;
    return exp(t)+sin(1.0);
}

double BorderParabolic1D2::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    double sgm = 3.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;
    sum += v1(t) * gause_a * exp(-((x-e1)*(x-e1))/gause_b);
    sum += v2(t) * gause_a * exp(-((x-e2)*(x-e2))/gause_b);
    sum += v3(t) * gause_a * exp(-((x-e3)*(x-e3))/gause_b);
    return sum;
}

double BorderParabolic1D2::v1(double t) const
{
    return 0.4*t;
}

double BorderParabolic1D2::v2(double t) const
{
    return 0.5*t;
}

double BorderParabolic1D2::v3(double t) const
{
    return 0.8*t;
}

