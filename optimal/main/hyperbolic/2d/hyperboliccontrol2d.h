#ifndef HYPERBOLICCONTROL2D_H
#define HYPERBOLICCONTROL2D_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>

class HyperbolicControl2D : public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter
{
public:
    HyperbolicControl2D();
    virtual ~HyperbolicControl2D();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double bfi1(unsigned int i) const;
    virtual double bfi2(unsigned int i) const;
    virtual double bm1(unsigned int j) const;
    virtual double bm2(unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    double v1(double t) const;
    double v2(double t) const;
    double fxt(double x, double t) const;

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    unsigned int N;
    unsigned int M;
    unsigned int L;
    double a;
    double alpha0;
    double alpha1;

    DoubleVector U0;
    DoubleVector U1;
};

#endif // HYPERBOLICCONTROL2D_H
