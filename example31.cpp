#include "example31.h"

// u(x,t)=t(x+1)

void BorderParabolic1D31::main(int argc, char* argv[])
{
    DoubleMatrix m;
    BorderParabolic1D31 bp;
    bp.calculateU(m, bp.hx, bp.ht, bp.N, bp.M, 2.0);
    IPrinter::printMatrix(m);
//    FILE *file = fopen("example311.txt", "w");
//    IPrinter::printVector(m[m.size()-1], NULL, bp.N, 0, 0, file);
//    fclose(file);
}

BorderParabolic1D31::BorderParabolic1D31()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);
    e1 = 0.30;
    e2 = 0.75;
}

double BorderParabolic1D31::initial(unsigned int i) const
{
    return 0.0;
}

double BorderParabolic1D31::boundary(Boundary type, unsigned int j) const
{
    double t = j*ht;
    if (type == Left)  return t;
    if (type == Right) return 2.0*t;
    return 0.0;
}

double BorderParabolic1D31::f(unsigned int i, unsigned int j) const
{
    double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    double sgm = 10.0*hx;
    double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    double gause_b = 2.0*sgm*sgm;
    sum += v1(t) * gause_a * exp(-((x-e1)*(x-e1))/gause_b);
    sum += v2(t) * gause_a * exp(-((x-e2)*(x-e2))/gause_b);
    return sum;
}

double BorderParabolic1D31::v1(double t) const
{
    return 20.0*t;
}

double BorderParabolic1D31::v2(double t) const
{
    return 30.0*t;
}
