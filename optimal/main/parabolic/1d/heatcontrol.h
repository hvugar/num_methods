#ifndef HEATCONTROL_H
#define HEATCONTROL_H

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

struct HeatControl : public RnFunction, public IGradient, public IParabolicEquation, public IBackwardParabolicEquation, public IPrinter
{
public:
    HeatControl();
    virtual ~HeatControl() {}

    virtual double fx(const DoubleVector& f);
    virtual void gradient(const DoubleVector& f, DoubleVector &g);

    inline virtual double fi(unsigned int i) const;
    inline virtual double m1(unsigned int j) const;
    inline virtual double m2(unsigned int j) const;
    inline virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi(unsigned int i) const;
    virtual double bm1(unsigned int j) const;
    virtual double bm2(unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;

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

    inline double u(double x, double t) const { return x*x+t*t; }
    inline double fxt(double x, double t) { return 2.0*t - 2.0*a; }
    const DoubleVector* pf;
    const DoubleVector* pu;
    DoubleVector U;

public:
    static void main();
};

#endif // HEATCONTROL_H
