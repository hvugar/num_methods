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

    e11 = 0.20;
    e12 = 0.30;

    e21 = 0.40;
    e22 = 0.70;

    e31 = 0.80;
    e32 = 0.50;
}

double BorderParabolic2D1::fi(unsigned int i, unsigned int j) const
{
    return 2.0;
}

double BorderParabolic2D1::m1(unsigned int j, unsigned int k) const
{
    double t = 0.5*k*ht;
    double x2 = j*h2;
    return 2.0 + t + 2.0*x2;
}

double BorderParabolic2D1::m2(unsigned int j, unsigned int k) const
{
    double t = 0.5*k*ht;
    double x2 = j*h2;
    return 2.0 + 0.8*t + 0.2*x2;
}

double BorderParabolic2D1::m3(unsigned int i, unsigned int k) const
{
    double t = 0.5*k*ht;
    double x1 = i*h1;
    return 2.0 + 1.2*t + sin(x1);
}

double BorderParabolic2D1::m4(unsigned int i, unsigned int k) const
{
    double t = 0.5*k*ht;
    double x1 = i*h1;
    return 2.0 + 1.5*t + 0.5*x1;
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

    sum += v1(t) * gause_a * exp(-((x1-e11)*(x1-e11) + (x2-e12)*(x2-e12))/gause_b);
    sum += v2(t) * gause_a * exp(-((x1-e21)*(x1-e21) + (x2-e22)*(x2-e22))/gause_b);
    sum += v3(t) * gause_a * exp(-((x1-e31)*(x1-e31) + (x2-e32)*(x2-e32))/gause_b);
    return sum;
}

double BorderParabolic2D1::v1(double t) const
{
    return 4.0*t;
}

double BorderParabolic2D1::v2(double t) const
{
    return 5.0*t;
}

double BorderParabolic2D1::v3(double t) const
{
    return 6.0*t;
}
