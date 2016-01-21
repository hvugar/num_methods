#ifndef POINTCONTROL1_H
#define POINTCONTROL1_H

#include <math.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <function.h>
#include <printer.h>

class PointControl1 : public RnFunction, public IPrinter
{
public:
    PointControl1(double t0, double t1, double x0, double x1, double dt, double dx);
    virtual ~PointControl1() {}

    virtual double fx(const DoubleVector &p);
    virtual void gradient(const DoubleVector& p, DoubleVector& g, double gradient_step);

    virtual void print(unsigned int iterationCount, const DoubleVector& x, const DoubleVector &g, double m_alpha, RnFunction* f) const;

    void calculate_x(const DoubleVector& p);
    void calculate_psi();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double dt;
    double dx;
    unsigned int n;
    DoubleVector T;

    double f(double t, double x);
    double dxdt(double t, double x, const DoubleVector& p);
    double px(double t, double psi, double x);
    double delta(double t);

    DoubleVector x;
    DoubleVector psi;
    DoubleVector p;

    void write(DoubleVector &x);

public:
    static void main();
};

#endif // POINTCONTROL1_H
