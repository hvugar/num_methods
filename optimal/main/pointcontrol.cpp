#include "pointcontrol.h"
#include <math.h>
#include "utils.h"

PointControl::PointControl(double t0, double t1, double x0, double x1, double dt, double dx)
{
    this->t0 = 0.0;
    this->t1 = 1.0;
    this->x0 = 0.0;
    this->x1 = 1.7391017563;
    this->dx = 0.00001;
    this->dt = 0.00001;

    n = (unsigned int)(ceil((t1 - t0)/dt))+1;
    epsilon = 0.02;

    T.resize(2);
    T[0] = 0.2;
    T[1] = 0.5;

    p.resize(2);
    p[0] = 0.3;
    p[1] = 0.4;
}

double PointControl::fx(const DoubleVector &p)
{
    DoubleVector x(n, 0.0);
    calculate_x(x);

    printX("x:",x);
    return 0.0;
}

void PointControl::gradient(double step, const DoubleVector &x, DoubleVector &g)
{
}

void PointControl::calculate()
{
}

double PointControl::f(double x, double t)
{
    return 3.0*t*t - t*t*t + x;
}

double PointControl::delta(double t)
{
//    for (unsigned int i=0; i<T.size(); i++)
//    {
//        if (fabs(T[i] - t) < ((epsilon / 2.0) + 0.000001)) return 1.0 / epsilon;
//    }
//    return 0.0;
    double dlt = 0.0;
    if (fabs(t - T[0]) <= ((epsilon / 2.0) + 0.000001)) dlt = 1.0 / epsilon;
    if (fabs(t - T[1]) <= ((epsilon / 2.0) + 0.000001)) dlt = 1.0 / epsilon;
    return dlt;
}

double PointControl::dxdt(double x, double t)
{
    double sum = f(x, t);
//    for (unsigned int i=0; i<p.size(); i++)
//    {
//        sum += p[i]*delta(t);
//    }
    return sum;
}

void PointControl::calculate_x(DoubleVector &x)
{
    int i=0;
    x[i] = x0;

    double _t0 = t0;
    double _t1 = t1;
    double _x0 = x0;
//    double _x1 = x1;

    while (_t0 <= _t1)
    {
        double k1 = dxdt(_t0, _x0);
        double k2 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k1);
        double k3 = dxdt(_t0+dt/2.0, _x0+(dt/2.0)*k2);
        double k4 = dxdt(_t0+dt, _x0+dt*k3);

        _x0 = _x0 + (dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        _t0 = _t0 + dt;

        i++;
        x[i] = _x0;
    }
}
