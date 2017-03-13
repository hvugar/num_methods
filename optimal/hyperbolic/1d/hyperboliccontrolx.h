#ifndef HYPERBOLICCONTROLX_H
#define HYPERBOLICCONTROLX_H

#include <stdio.h>
#include <math.h>
#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>
#include <gradient_cjt.h>

class HyperbolicControlX : public R1Function, public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter
{
public:
    HyperbolicControlX();
    virtual ~HyperbolicControlX() {}

    virtual double fx(double t) const;
    virtual double fx(const DoubleVector& v) const;
    virtual void gradient(const DoubleVector &v, DoubleVector &g);

    virtual double initial1(unsigned int i) const;
    virtual double initial2(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual double binitial1(unsigned int i) const;
    virtual double binitial2(unsigned int i) const;
    virtual double bboundary(Boundary type, unsigned int j) const;
    virtual double bf(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int i, const DoubleVector& f0, const DoubleVector &g, double fx, GradientMethod::MethodResult result) const;

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

#endif // HYPERBOLICCONTROLX_H
