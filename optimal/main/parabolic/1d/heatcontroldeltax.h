#ifndef HEATCONTROLDELTAX_H
#define HEATCONTROLDELTAX_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

/**
 * @brief The HeatControl struct
 * du/dt = a(d^2u/dx^2) + f(x,t);
 * u(x,0) = fi(x);
 * u(0,t) = m1(t);
 * u(l,t) = m2(t);
 */

struct HeatControlDeltaX : public RnFunction, public IGradient, public IParabolicEquation, public IBackwardParabolicEquation, public IPrinter, public Projection
{
public:
    HeatControlDeltaX();
    virtual ~HeatControlDeltaX() {}

    virtual double fx(const DoubleVector& e);
    virtual void gradient(const DoubleVector& e, DoubleVector &g);

    inline virtual double fi(unsigned int i) const;
    inline virtual double m1(unsigned int j) const;
    inline virtual double m2(unsigned int j) const;
    inline virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi(unsigned int i) const;
    virtual double bm1(unsigned int j) const;
    virtual double bm2(unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
    virtual void project(DoubleVector &x, int index);
private:
    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    unsigned int N;
    unsigned int M;
    double a;
    unsigned int L;

    double f1(double t) const { return 2.0*t; }
    inline double u(double x, double t) const { return x*x+t*t; }
    const DoubleVector* pe;
    const DoubleVector* pu;
    DoubleVector U;

public:
    static void main();
};

#endif // HEATCONTROLDELTAX_H
