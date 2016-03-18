#ifndef PARABOLICCONTROL1D332_H
#define PARABOLICCONTROL1D332_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <gradient.h>
#include <function.h>
#include <gradient_cjt.h>
#include <math.h>

// u(x,t) = 0.5*x + t^2

class Parabolic1DControl332 : public IParabolicEquation, public IBackwardParabolicEquation, public RnFunction, public IGradient, public IPrinter, public Projection
{
public:
    Parabolic1DControl332();
    virtual ~Parabolic1DControl332() {}

    virtual double fx(const DoubleVector &e);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &v, const DoubleVector &gradient, double alpha, RnFunction *fn) const;
    virtual void project(DoubleVector &x, int index);

    virtual double fi(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi(unsigned int i) const;
    virtual double bm1(unsigned int j) const;
    virtual double bm2(unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    static void main();

private:
    double x0;
    double x1;
    double t0;
    double t1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    unsigned int L;
    double a;

    double v1(double t) const;
    double v2(double t) const;

    DoubleVector U;
    const DoubleVector *pe;
    const DoubleVector *pu;
};

#endif // PARABOLICCONTROL1D332_H
