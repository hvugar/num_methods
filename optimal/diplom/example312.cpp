#include "example312.h"

// u(x,t)=sin(x)+e^t

void BorderParabolic1D312::main()
{
    DoubleMatrix m;
    BorderParabolic1D312 bp;
    bp.calculateU(m, bp.hx, bp.ht, bp.N, bp.M, 2.0);
    IPrinter::printMatrix(m);
    FILE *file = fopen("example312.txt", "w");
    IPrinter::printVector(m[m.size()-1], NULL, bp.N, 0, 0, file);
    fclose(file);
}

BorderParabolic1D312::BorderParabolic1D312()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);
    L = 3;
    e.resize(L);
    e[0] = 0.20;
    e[1] = 0.55;
    e[2] = 0.84;
}

double BorderParabolic1D312::fi(unsigned int i) const
{
    double x = i*hx;
    return sin(x)+1.0;
}

double BorderParabolic1D312::m1(unsigned int j) const
{
    double t = j*ht;
    return exp(t);
}

double BorderParabolic1D312::m2(unsigned int j) const
{
    double t = j*ht;
    return exp(t)+sin(1.0);
}

double BorderParabolic1D312::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    double sgm = 3.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;
    sum += v1(t) * gause_a * exp(-((x-e[0])*(x-e[0]))/gause_b);
    sum += v2(t) * gause_a * exp(-((x-e[1])*(x-e[1]))/gause_b);
    sum += v3(t) * gause_a * exp(-((x-e[2])*(x-e[2]))/gause_b);
    return sum;
}

double BorderParabolic1D312::v1(double t) const
{
    return 50.0*t;
}

double BorderParabolic1D312::v2(double t) const
{
    return 10.0*t;
}

double BorderParabolic1D312::v3(double t) const
{
    return 40.0*t;
}

