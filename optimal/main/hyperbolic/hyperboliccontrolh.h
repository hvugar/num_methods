#ifndef HYPERBOLICCONTROLH_H
#define HYPERBOLICCONTROLH_H

#include <stdio.h>
#include <math.h>
#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <tomasmethod.h>

class HyperbolicControlH : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter
{
public:
    HyperbolicControlH();
    virtual ~HyperbolicControlH() {}

    virtual double fx(double t);
    virtual double fx(const DoubleVector& v);
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

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

    virtual void print(unsigned int iteration, const DoubleVector& v, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    static void main();

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double hx;
    double ht;
    double a;
    unsigned int N;
    unsigned int M;
    unsigned int D;
    unsigned int L;
    unsigned int Xi;
    double xi;
    const DoubleVector *pv;
    const DoubleMatrix *pu;
    double U;
    double lamda;
};

#endif // HYPERBOLICCONTROLH_H
