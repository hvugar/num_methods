#include "cfunction4.h"

void ControlFunction4::Main(int argc, char *agrv[])
{
    ControlFunction4 c;
    DoubleVector t(c.N+1);
    DoubleVector x1(c.N+1);
    DoubleVector x2(c.N+1);
    DoubleVector u1(c.N+1);
    DoubleVector u2(c.N+1);

    for (unsigned int i=0; i<=c.N; i++)
    {
        t[i] = i*c.h;
        u1[i] = 2.0*t[i];
        u2[i] = 3.0*t[i];
        x1[i] = x2[i] = 0.0;
    }
    c.calculate_X(u1, u2, x1, x2);

    IPrinter::printVector(t, "t ");
    IPrinter::printVector(x1, "x1");
    IPrinter::printVector(x2, "x2");
    IPrinter::printVector(u1, "u1");
    IPrinter::printVector(u2, "u2");
}

ControlFunction4::ControlFunction4()
{
    t0 = 0.0;
    t1 = 1.0;
    N = 1000;
    h = (t1 - t0) / N;
}

ControlFunction4::~ControlFunction4()
{}

double ControlFunction4::fx(const DoubleVector &x)
{
    return 0.0;
}

void ControlFunction4::gradient(const DoubleVector &x, DoubleVector &g)
{}

void ControlFunction4::print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const
{}

void ControlFunction4::calculate_X(const DoubleVector &u1, const DoubleVector &u2, DoubleVector &x1, DoubleVector &x2)
{
    double k1[] = {0.0, 0.0};
    double k2[] = {0.0, 0.0};
    double k3[] = {0.0, 0.0};
    double k4[] = {0.0, 0.0};

    x1[0] = 0.0;
    x2[0] = 0.0;
//    DoubleVector _x(2);
//    _x[0] = x1[0];
//    _x[1] = x2[0];

    h = +fabs(h);

    for (unsigned int i=0; i<N; i++)
    {
//        _x[0] = x1[i];
//        _x[1] = x2[i];

        double t = i*h;

        k1[0] = fx1(t, x1[i], x2[i], u1[i], u2[i]);
        k1[1] = fx2(t, x1[i], x2[i], u1[i], u2[i]);

//        _x[0] = x1[i] + (h/2.0) * k1[0];
//        _x[1] = x2[i] + (h/2.0) * k1[1];
        k2[0] = fx1(t+h/2.0, x1[i] + (h/2.0) * k1[0], x2[i] + (h/2.0) * k1[1], u1[i], u2[i]);
        k2[1] = fx2(t+h/2.0, x1[i] + (h/2.0) * k1[0], x2[i] + (h/2.0) * k1[1], u1[i], u2[i]);

//        _x[0] = x1[i] + (h/2.0) * k2[0];
//        _x[1] = x2[i] + (h/2.0) * k2[1];
        k3[0] = fx1(t+h/2.0, x1[i] + (h/2.0) * k2[0], x2[i] + (h/2.0) * k2[1], u1[i], u2[i]);
        k3[1] = fx2(t+h/2.0, x1[i] + (h/2.0) * k2[0], x2[i] + (h/2.0) * k2[1], u1[i], u2[i]);

//        _x[0] = x1[i] + h * k3[0];
//        _x[1] = x2[i] + h * k3[1];
        k4[0] = fx1(t+h, x1[i] + h * k3[0], x2[i] + h * k3[1], u1[i], u2[i]);
        k4[1] = fx2(t+h, x1[i] + h * k3[0], x2[i] + h * k3[1], u1[i], u2[i]);

        x1[i+1] = x1[i] + (h/6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
        x2[i+1] = x2[i] + (h/6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
    }
}

double ControlFunction4::fx1(double t, double x1, double x2, double u1, double u2)
{
    return 3.0*x1 + 2.0*x2 + u1 - (3.0*t*t + 2.0*t*t*t);
}

double ControlFunction4::fx2(double t, double x1, double x2, double u1, double u2)
{
    return 4.0*x1 + t*x2 + 2.0*u1 + 3.0*u2 - (t*t + t*t*t*t + 13.0*t);
}
