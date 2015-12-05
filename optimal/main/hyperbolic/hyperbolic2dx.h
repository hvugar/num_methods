#ifndef HYPERBOLIC2DX_H
#define HYPERBOLIC2DX_H

#include <function.h>
#include <printer.h>
#include <projection.h>

struct Hyperbolic2DX : public RnFunction, Projection, Printer
{
    Hyperbolic2DX();
    virtual ~Hyperbolic2DX() {}

    virtual double fx(const DoubleVector& x);
    virtual void gradient(const DoubleVector& x, DoubleVector& g, double gradient_step=0.000001);

    virtual void project(DoubleVector &x, int index);

    virtual void print(unsigned int iteration, const DoubleVector& x, const DoubleVector &gradient, double alpha, RnFunction* fn);
};

#endif // HYPERBOLIC2DX_H
