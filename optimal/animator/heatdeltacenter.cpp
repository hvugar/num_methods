#include "heatdeltacenter.h"

void HeatDeltaCenter::main()
{
    HeatDeltaCenter hc;

    DoubleMatrix u;
    hc.calculateU(u, hc.hx, hc.hy, hc.ht, hc.N1, hc.N2, hc.M);
}

HeatDeltaCenter::HeatDeltaCenter()
{
    t0 = 0.0;
    t1 = 1.0;
    x0 = 0.0;
    x1 = 1.0;
    y0 = 0.0;
    y1 = 1.0;
    M = 1000;
    N1 = 1000;
    N2 = 1000;
    ht = 0.001;
    hx = 0.001;
    hy = 0.001;
}

double HeatDeltaCenter::fi(unsigned int i, unsigned int j) const
{
    return 0.0;
}

double HeatDeltaCenter::m1(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HeatDeltaCenter::m2(unsigned int j, unsigned int k) const
{
    return 0.0;
}

double HeatDeltaCenter::m3(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HeatDeltaCenter::m4(unsigned int i, unsigned int k) const
{
    return 0.0;
}

double HeatDeltaCenter::f(unsigned int i, unsigned int j, unsigned int k) const
{
    //    unsigned int k1= k%2==0 ? k/2 : (k+1)/2;

    double x1 = i*hx;
    double x2 = j*hy;
    double t  = 0.5*k*ht;
    double sum = 0.0;

    double sgm1 = 3.0*hx;
    double sgm2 = 3.0*hy;
    double gause_a = 1.0/(2.0*M_PI*sgm1*sgm2);
    double gause_b = 2.0*sgm1*sgm2;

    sum += v1(t) * gause_a * exp(-((x1-0.5)*(x1-0.5) + (x2-0.5)*(x2-0.5))/gause_b);
    return sum;
}

double HeatDeltaCenter::v1(double t) const
{
    return 10.0*t;
}

