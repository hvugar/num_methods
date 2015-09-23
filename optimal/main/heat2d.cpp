#include "heat2d.h"

Heat2DControl::Heat2DControl()
{
    N = 1000;
    M = 1000;
    t0 = 0.0;
    t1 = 1.0;

    x11 = 0.0;
    x12 = 1.0;
    x21 = 0.0;
    x22 = 1.0;
}

Heat2DControl::~Heat2DControl()
{

}

double Heat2DControl::fx(const DoubleVector &x)
{
    return 0.0;
}

void Heat2DControl::gradient(double step, const DoubleVector &x, DoubleVector &g)
{

}

double Heat2DControl::f(double x1, double x2, double t) { return 2.0*t - 4.0; }
double Heat2DControl::fi(double x1, double x2) { return x1*x1 + x2*x2; }
double Heat2DControl::m1(double x2, double t) { return x2*x2 + t*t; }
double Heat2DControl::m2(double x2, double t) { return x2*x2 + t*t + 1.0; }
double Heat2DControl::m3(double x1, double t) { return x1*x1 + t*t; }
double Heat2DControl::m4(double x1, double t) { return x1*x1 + t*t + 1.0; }

