#include "example1.h"

double OrdDifEquationA::fx(double t, double a) const
{
    return t*a + (2.0*t-t*t*t-0.75);
}

double OrdDifEquationB::fx(double t, double b) const
{
    return t*b + 3.0;
}

double OrdDifEquationX::fx(double t, double x) const
{
    double K = 0.0;
    return t*x - 3.0*K*x;
}

Example1::Example1()
{
    t0 = 0.0;
    t1 = 1.0;
    h = 0.001;
    N = 1000;

//    double b = 0.0;

    alphaOde ao;
    ao.e = this;
    bettaOde bo;
    bo.e = this;

    //////////////////////////////////////

    DoubleVector a;
    DoubleVector b;

    runge_kutta(&ao, 2.0, 0.0, N, a, h);
    runge_kutta(&bo, 0.0, 0.0, N, b, h);

    IPrinter::printVector(a, "alpha");
    IPrinter::printVector(b, "betta");

    DoubleVector x(N+1);

    x[500] = a[500] / (1.0 - b[500]);

    for (unsigned int i=0; i<=N; i++)
    {
        x[i] = a[i] + b[i]*x[500];
    }
    IPrinter::printVector(x, "x    ");
}

void Example1::runge_kutta(OrdDifEquation *ode, double x0, double y0, unsigned int N, DoubleVector &y, double h)
{
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double k4 = 0.0;

    y.clear();
    y.resize(N+1,0.0);

    y[0] = y0;
    for (unsigned int i=0; i<N; i++)
    {
        double _x = x0 + i*h;
        double _y = y[i];

        k1 = ode->fx(_x, _y);
        k2 = ode->fx(_x+h/2.0, _y+(h/2.0)*k1);
        k3 = ode->fx(_x+h/2.0, _y+(h/2.0)*k2);
        k4 = ode->fx(_x+h, _y+h*k3);
        y[i+1] = y[i] + (h/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);;
    }
}

double Example1::getK(double x1, double x2)
{
    if (x1 >= DBL_MIN && x2 < -3.0) return 1.2;
    if (x1 >= -3.0    && x2 < +1.0) return 2.0;
    if (x1 >= +1.0    && x2 < +2.5) return 0.5;
    if (x1 >= +2.5    && x2 < +4.2) return 3.2;
    if (x1 >= +4.2    && x2 < +5.3) return 1.5;
    if (x1 >= +5.3 && x2 < DBL_MAX) return 0.8;
    return 0.0;
}

double Example1::A(double t) const
{
    return t;
}

double Example1::B(double t UNUSED_PARAM, unsigned int k) const
{
    double K = 0.5;
    if (k==1) return -3.0*K;
    return 0.0;
}

double Example1::C(double t UNUSED_PARAM) const
{
    return 2.0;
}
