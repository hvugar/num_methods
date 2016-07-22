#ifndef EXAMPLE1_H
#define EXAMPLE1_H

#include <function.h>
#include <printer.h>
#include <stdlib.h>
#include <float.h>


class ODELoadedSystem
{
public:
    unsigned int S;

    virtual double A(double t) const = 0;
    virtual double B(double t, unsigned int k) const = 0;
    virtual double C(double t) const = 0;
};

struct OrdDifEquationA : public OrdDifEquation
{
    virtual double fx(double t, double a) const;
};

struct OrdDifEquationB : public OrdDifEquation
{
    virtual double fx(double t, double b) const;
};

class Example1 : public ODELoadedSystem
{
public:
    Example1();

    void runge_kutta(OrdDifEquation *ode, double x0, double y0, unsigned int N, DoubleVector &y, double h);

    double t0;
    double t1;
    double h;
    unsigned int N;

    double getK(double x1, double x2);

    virtual double A(double t) const;
    virtual double B(double t, unsigned int k) const;
    virtual double C(double t) const;

    struct alphaOde : public OrdDifEquation {
        Example1 *e;
        virtual double fx(double t, double a) const {
            return e->A(t)*a + e->C(t);
        }
    };

    struct bettaOde : public OrdDifEquation {
        Example1 *e;
        virtual double fx(double t, double b) const {
            return e->A(t)*b + e->B(t, 1);
        }
    };
};

struct OrdDifEquationX : public OrdDifEquation
{
    virtual double fx(double t, double x) const;
    Example1 *e;
};

#endif // EXAMPLE1_H
