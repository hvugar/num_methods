#include "borderhyperbolic2d.h"

#define SAMPLE6

double k1;
double k2;
double k3;
double k4;
double k5;
double k6;

BorderHyperbolic2D::BorderHyperbolic2D()
{
    a1 = a2 = 1.0;
    x10 = x20 = t0 = 0.0;
    x11 = x21 = 1.0;
    t1 = 20.0;
    h1 = 0.010;
    h2 = 0.010;
    ht = 0.005;
    N1 = (unsigned int)(ceil(x11-x10)/h1);
    N2 = (unsigned int)(ceil(x21-x20)/h2);
    M  = (unsigned int)(ceil(t1-t0)/ht);

    k1 = 50.0;
    k2 = 10.0;
    k3 = 30.0;
    k4 = 50.0;
    k5 = 40.0;
    k6 = 50.0;

    qamma = 0.0;
}

BorderHyperbolic2D::~BorderHyperbolic2D()
{
}

double BorderHyperbolic2D::u(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
#ifdef SAMPLE1
    return x1*x1 + x2*x2 + t*t;
#endif
#ifdef SAMPLE2
    return x1*x1 + x1*x2 + sin(x1) + cos(x2) + exp(t);
#endif
#ifdef SAMPLE3
    return sin(x1)*sin(x1) + cos(x2) + x1*x2 + t*t*t;
#endif
#ifdef SAMPLE4
    return sin(t) + exp(x1) + x2*x2;
#endif
#ifdef SAMPLE5
    //return x1*x1 + x1*x2 + 5.0*sin(10.0*x1) + 3.0*cos(5.0*x2) + 4.0*cos(5.0*x1*t);
    return x1*x1 + x1*x2 + k1*sin(k2*x1) + k3*cos(k4*x2) + k5*cos(k6*x1*t);
#endif
#ifdef SAMPLE6
    return 2.0*sin(3.0*x1) + 4.0*cos(2.0*x2) + sin(2.0*t);
#endif
}

double BorderHyperbolic2D::fi1(unsigned int i, unsigned int j) const
{
    return u(i, j, 0);
}

double BorderHyperbolic2D::fi2(unsigned int i, unsigned int j) const
{
#ifdef SAMPLE1
    return 0.0;
#endif
#ifdef SAMPLE2
    return 1.0;
#endif
#ifdef SAMPLE3
    return 0.0;
#endif
#ifdef SAMPLE4
    return 1.0;
#endif
#ifdef SAMPLE5
    return 0.0;
#endif
#ifdef SAMPLE6
    return 2.0;
#endif
}

double BorderHyperbolic2D::m1(unsigned int j, unsigned int k) const
{
    return u(0, j, k);
}

double BorderHyperbolic2D::m2(unsigned int j, unsigned int k) const
{
    return u(N1, j, k);
}

double BorderHyperbolic2D::m3(unsigned int i, unsigned int k) const
{
    return u(i, 0, k);
}

double BorderHyperbolic2D::m4(unsigned int i, unsigned int k) const
{
    return u(i, N2, k);
}
double BorderHyperbolic2D::f(unsigned int i, unsigned int j, unsigned int k) const
{
    double x1 = i*h1;
    double x2 = j*h2;
    double t  = k*ht;
    C_UNUSED(x1);
    C_UNUSED(x2);
    C_UNUSED(t);
#ifdef SAMPLE1
    return 2.0 - 2.0*(a1*a1) - 2.0*(a2*a2) + qamma*2.0*t;
#endif
#ifdef SAMPLE2
    return exp(t) + (a1*a1)*(sin(x1) - 2.0) + (a2*a2)*cos(x2);// + qamma*exp(t);
#endif
#ifdef SAMPLE3
    return 6*t - (a1*a1)*(2.0*cos(2.0*x1)) + (a2*a2)*cos(x2) ;
#endif
#ifdef SAMPLE4
    return -sin(t) - (a1*a1)*exp(x1) - (a2*a2)*2.0;
#endif
#ifdef SAMPLE5
    //return -100.0*x1*x1*cos(5.0*x1*t) + a1*a1*(500.0*sin(10.0*x1) + 100.0*t*t*cos(5.0*x1*t) - 2.0) + a2*a2*(75.0*cos(5.0*x2));
    return -k5*k6*k6*x1*x1*cos(k6*x1*t)
            - a1*a1*(2.0-k1*k2*k2*sin(k2*x1) - k5*k6*k6*t*t*cos(k6*x1*t))
            - a2*a2*(-k3*k4*k4*cos(k4*x2));
#endif
#ifdef SAMPLE6
    return -4.0*sin(2.0*t) + a1*a1*(18.0*sin(3.0*x1)) + a2*a2*(16.0*cos(2.0*x2)) + qamma*2.0*cos(2.0*t);
#endif
}


