#include "borderparabolic1d.h"

#define EXAMPLE2

void BorderParabolic1D::main()
{
    DoubleMatrix m;
    BorderParabolic1D bp;
    bp.calculateU(m, bp.hx, bp.ht, bp.N, bp.M, 2.0);
    IPrinter::printMatrix(m);
#ifdef EXAMPLE1
    FILE *file = fopen("example1.txt", "w");
#endif
#ifdef EXAMPLE2
    FILE *file = fopen("example2.txt", "w");
#endif
    //IPrinter::printMatrix(m, bp.N, bp.M, NULL, file);
    IPrinter::printVector(m[m.size()-1], NULL, bp.N, 0, 0, file);
    fclose(file);
}

BorderParabolic1D::BorderParabolic1D()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);

#ifdef EXAMPLE1
    e1 = 0.30;
    e2 = 0.75;
#endif
#ifdef EXAMPLE2
    e1 = 0.20;
    e2 = 0.55;
    e3 = 0.84;
#endif
}

BorderParabolic1D::~BorderParabolic1D()
{}

double BorderParabolic1D::fi(unsigned int i) const
{
    double x = i*hx;
#ifdef EXAMPLE1
    return 0.0;
#endif
#ifdef EXAMPLE2
    return sin(x)+1.0;
#endif
}

double BorderParabolic1D::m1(unsigned int j) const
{
    double t = j*ht;
#ifdef EXAMPLE1
    return 2.0*t;
#endif
#ifdef EXAMPLE2
    return exp(t);
#endif
}

double BorderParabolic1D::m2(unsigned int j) const
{
    double t = j*ht;
#ifdef EXAMPLE1
    return 1.0*t;
#endif
#ifdef EXAMPLE2
    return exp(t)+sin(1.0);
#endif
}

double BorderParabolic1D::f(unsigned int i, unsigned int j) const
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

double BorderParabolic1D::v1(double t) const
{
#ifdef EXAMPLE1
    return 2.0*t;
#endif
#ifdef EXAMPLE2
    return 0.4*t;
#endif
}

double BorderParabolic1D::v2(double t) const
{
#ifdef EXAMPLE1
    return 3.0*t;
#endif
#ifdef EXAMPLE2
    return 0.5*t;
#endif
}

double BorderParabolic1D::v3(double t) const
{
#ifdef EXAMPLE2
    return 0.8*t;
#endif
}

