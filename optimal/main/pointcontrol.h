#ifndef POINTCONTROL_H
#define POINTCONTROL_H

#include <function.h>

class PointControl : public RnFunction
{
public:
    PointControl(double t0, double t1, double x0, double x1, double dt, double dx);

    virtual double fx(const DoubleVector& p);
    virtual void gradient(double step, const DoubleVector& p, DoubleVector& g);

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

    double f(double t, double x);
    double delta(double t);
    double dxdt(double t, double x);
    double px(double t, double psi, double x);

    void calculate_x(DoubleVector &x, const DoubleVector& p);
    void calculate_psi(DoubleVector &x, DoubleVector &p);

public:
    static void main();
};

#endif // POINTCONTROL_H
