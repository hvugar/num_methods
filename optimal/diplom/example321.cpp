#include "example321.h"

// u(x,t)=sin(x)+e^t

void Parabolic1DControl321::main()
{
//    DoubleMatrix m;
    Parabolic1DControl321 bp;
//    bp.calculateU(m, bp.hx, bp.ht, bp.N, bp.M, 2.0);
//    IPrinter::printMatrix(m);
//    FILE *file = fopen("example321.txt", "w");
//    IPrinter::printVector(m[m.size()-1], NULL, bp.N, 0, 0, file);
//    fclose(file);
}

Parabolic1DControl321::Parabolic1DControl321()
{
    x0 = 0.0;
    x1 = 1.0;
    t0 = 0.0;
    t1 = 1.0;
    hx = 0.001;
    ht = 0.001;
    N = (unsigned int)(ceil(x1-x0)/hx);
    M = (unsigned int)(ceil(t1-t0)/ht);
    a = 1.0;
    e.resize(L);
    e[0] = 0.10;
    e[1] = 0.85;

    IParabolicEquation::calculateU(U, hx, ht, N, M, a);
    IPrinter::printVector(U);
    FILE *file = fopen("example321.txt", "w");
    IPrinter::printVector(U, NULL, N, 0, 0, file);
    fclose(file);
}

double Parabolic1DControl321::fx(const DoubleVector &x)
{
    return 0.0;
}

double Parabolic1DControl321::fi(unsigned int i) const
{
    double x = i*hx;
    return x+cos(x);
}

double Parabolic1DControl321::m1(unsigned int j) const
{
    double t = j*ht;
    return t*t+1.0;
}

double Parabolic1DControl321::m2(unsigned int j) const
{
    double t = j*ht;
    return t*t+1.0+cos(1.0);
}

double Parabolic1DControl321::f(unsigned int i, unsigned int j) const
{
    //double x = i*hx;
    double t = j*ht;
    double sum = 0.0;
    //double sgm = 3.0*hx;
    //double gause_a = 1.0/(sqrt(2.0*M_PI)*sgm);
    //double gause_b = 2.0*sgm*sgm;
    //sum += v1(t) * gause_a * exp(-((x-e[0])*(x-e[0]))/gause_b);
    //sum += v2(t) * gause_a * exp(-((x-e[1])*(x-e[1]))/gause_b);

    if (i==e[0]/hx) sum += (1/hx)*v1(t);
    if (i==e[1]/hx) sum += (1/hx)*v2(t);

    return sum;
}

double Parabolic1DControl321::v1(double t) const
{
    return 102.0*t;
}

double Parabolic1DControl321::v2(double t) const
{
    return 100.0*t;
}

