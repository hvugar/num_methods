#include "example1.h"

double OrdDifEquationA::fx(double t, double a) const
{
    return t*a + (2.0*t-t*t*t-0.75);
}

double OrdDifEquationB::fx(double t, double b) const
{
    return t*b + 3.0;
}

Example1::Example1()
{
    double t0 = 0.0;
    unsigned int N = 1000;
    double h = 0.001;
    DoubleVector a;
    DoubleVector b;

    OrdDifEquationA odeA;
    OrdDifEquationB odeB;

    runge_kutta(&odeA, 0.0, 0.0, N, a, h);
    runge_kutta(&odeB, 0.0, 0.0, N, b, h);

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
