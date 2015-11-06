#ifndef HYPERBOLICCONTROL1D_H
#define HYPERBOLICCONTROL1D_H

#include <function.h>
#include <hyperbolicequation.h>
#include <doublevector.h>
#include <printer.h>
#include <projection.h>

class HyperbolicControl1D : public RnFunction, HyperbolicEquation, Printer
{
public:
    HyperbolicControl1D();
    virtual ~HyperbolicControl1D();

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual double fi1(unsigned int i) const;
    virtual double fi2(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn) const;

    static void main();

protected:
    double dt;
};

#endif
