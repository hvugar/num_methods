#ifndef HEATCONTROLDELTAF_H
#define HEATCONTROLDELTAF_H

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

struct HeatControlDeltaF : public RnFunction, public IGradient, public IParabolicEquation, public IBackwardParabolicEquation, public IPrinter
{
public:
    HeatControlDeltaF();
    virtual ~HeatControlDeltaF() {}

    virtual double fx(const DoubleVector& f) const;
    virtual void gradient(const DoubleVector& f, DoubleVector &g);

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double binitial(unsigned int i) const;
    virtual double bboundary(Boundary type, unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;
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
    double e;
    unsigned int E;

    double f1(double t) const { return 2.0*t+3.0; }
    inline double u(double x, double t) const { return x*x+t*t; }
    const DoubleVector* pf;
    const DoubleVector* pu;
    DoubleVector U;

public:
    static void main(int argc, char ** argv);
};

#endif // HEATCONTROLDELTAF_H
