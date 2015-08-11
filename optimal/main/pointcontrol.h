#ifndef POINTCONTROL_H
#define POINTCONTROL_H

#include <function.h>

class PointControl : public RnFunction
{
public:
    PointControl(double t0, double t1, double x0, double x1, double dt, double dx);

    virtual double fx(const DoubleVector& p);
    virtual void gradient(double step, const DoubleVector& x, DoubleVector& g);

    void calculate();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double dt;
    double dx;
    unsigned int n;
    double epsilon;

    DoubleVector T;
    DoubleVector p;

    double f(double x, double t);
    double delta(double t);
    double dxdt(double x, double t);

    void calculate_x(DoubleVector &x);
};

#endif // POINTCONTROL_H
