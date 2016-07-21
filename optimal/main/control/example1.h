#ifndef EXAMPLE1_H
#define EXAMPLE1_H

#include <function.h>
#include <printer.h>

struct OrdDifEquationA : public OrdDifEquation
{
    virtual double fx(double t, double a) const;
};


struct OrdDifEquationB : public OrdDifEquation
{
    virtual double fx(double t, double b) const;
};

class Example1
{
public:
    Example1();

    void runge_kutta(OrdDifEquation *ode, double x0, double y0, unsigned int N, DoubleVector &y, double h);
};

#endif // EXAMPLE1_H
