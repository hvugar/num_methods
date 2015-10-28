#ifndef HEATCONTROL1_H
#define HEATCONTROL1_H

#include <function.h>
#include <parabolicequation.h>

struct HeatControl1 : public RnFunction, ParabolicEquation
{
    HeatControl1();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fi(unsigned int i, double x) const;
    virtual double m1(unsigned int j, double t) const;
    virtual double m2(unsigned int j, double t) const;
    virtual double f(unsigned int i, double x, unsigned int j, double t) const;

    double u(double x, double t) const;

    void calculate(const DoubleVector& f);
};

#endif // HEATCONTROL1_H
