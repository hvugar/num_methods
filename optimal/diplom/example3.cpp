#include "example3.h"

void BorderParabolic2D1::main()
{
    DoubleMatrix m;
    BorderParabolic2D1 bp;
    bp.caluclateMVD(m, bp.h1, bp.h2, bp.ht, bp.N1, bp.N2, bp.M, 1.0, 1.0);
    IPrinter::printMatrix(m);
    FILE *file = fopen("example3.txt", "w");
    IPrinter::printMatrix(m, bp.N2, bp.N1, NULL, file);
    fclose(file);
}

BorderParabolic2D1::BorderParabolic2D1()
{
    x10 = 0.0;
    x11 = 1.0;

    x20 = 0.0;
    x21 = 1.0;

    t0 = 0.0;
    t1 = 1.0;

    h1 = 0.005;
    h2 = 0.005;
    ht = 0.005;

    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M = (unsigned int)(ceil(t1-t0)/ht);
    S = 3;

    e.resize(2*S);
    e[0] = 0.20;
    e[1] = 0.30;
    e[2] = 0.40;
    e[3] = 0.70;
    e[4] = 0.80;
    e[5] = 0.50;
}

double BorderParabolic2D1::fi(unsigned int i, unsigned int j) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    return x1*x2;
}

double BorderParabolic2D1::m1(unsigned int j, unsigned int k) const
{
    double x1 = 0.0*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    return x1*t + x2*t*t + x1*x2;
}

double BorderParabolic2D1::m2(unsigned int j, unsigned int k) const
{
    double x1 = N1*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    return x1*t + x2*t*t + x1*x2;
}

double BorderParabolic2D1::m3(unsigned int i, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = 0.0*h2;
    double t  = 0.5*k*ht;
    return x1*t + x2*t*t + x1*x2;
}

double BorderParabolic2D1::m4(unsigned int i, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = N2*h2;
    double t  = 0.5*k*ht;
    return x1*t + x2*t*t + x1*x2;
}

double BorderParabolic2D1::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = 0.5*k*ht;
    double sum = 0.0;

    double sgm1 = 3.0*h1;
    double sgm2 = 3.0*h2;
    double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    double gause_b = 2.0*sgm1*sgm2;

    sum += v1(t) * gause_a * exp(-((x1-e[0])*(x1-e[0]) + (x2-e[1])*(x2-e[1]))/gause_b);
    sum += v2(t) * gause_a * exp(-((x1-e[2])*(x1-e[2]) + (x2-e[3])*(x2-e[3]))/gause_b);
    sum += v3(t) * gause_a * exp(-((x1-e[4])*(x1-e[4]) + (x2-e[5])*(x2-e[5]))/gause_b);
    return sum;
}

double BorderParabolic2D1::v1(double t) const
{
    return 8.0*t;
}

double BorderParabolic2D1::v2(double t) const
{
    return 5.0*t;
}

double BorderParabolic2D1::v3(double t) const
{
    return 6.0*t;
}
