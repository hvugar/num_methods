#ifndef HEATCONTROLDELTAX_H
#define HEATCONTROLDELTAX_H

#include <math.h>
#include <stdlib.h>
#include <function.h>
#include <parabolicequation.h>
#include <printer.h>
#include <projection.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>

/**
 * @brief The HeatControl struct
 * du/dt = a(d^2u/dx^2) + f(x,t);
 * u(x,0) = fi(x);
 * u(0,t) = m1(t);
 * u(l,t) = m2(t);
 */

struct HeatControlDeltaX : public RnFunction, public IGradient, public IParabolicEquation,
        public IBackwardParabolicEquation, public IPrinter, public IProjection
{
public:
    HeatControlDeltaX();
    virtual ~HeatControlDeltaX() {}

    virtual double fx(const DoubleVector& e);
    virtual void gradient(const DoubleVector& e, DoubleVector &g);

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    inline virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi(unsigned int i) const;
    virtual double bm1(unsigned int j) const;
    virtual double bm2(unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &s, double a, RnFunction* f) const;
    virtual void project(DoubleVector &x, int index);
private:
    inline double u(double x, double t) const { return x*x+t*t; }
    double norm(const DoubleVector& v) const;

    double t0;
    double t1;
    double x0;
    double x1;
    double ht;
    double hx;
    double a;
    unsigned int N;
    unsigned int M;
    unsigned int L;

    double v1(double t) const { return 20.0*t; }
    double v2(double t) const { return 30.0*t; }
    double v3(double t) const { return 25.0*t; }

    const DoubleVector* pe;
    const DoubleVector* pu;
    DoubleVector U;

public:
    static void main();
};

#endif // HEATCONTROLDELTAX_H
