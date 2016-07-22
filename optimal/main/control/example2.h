#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <function.h>
#include <printer.h>
#include <stdlib.h>
#include <float.h>
#include <function.h>
#include <gradient.h>
#include <gradient_cjt.h>

class Example2  : public OrdDifEquation, public ConjugateGradient, public RnFunction, public IGradient, public IPrinter
{
public:
    Example2();

    virtual double fx(double t, double x) const;

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;

    static void main(int argc, char* argv[]);

    void calculateX(DoubleVector &x);

private:
    double t0;
    double t1;
    double h;
    unsigned int M;

    double x0;

    const DoubleVector *pK;
    const DoubleVector *px;

    double getK(double x) const;
    void runge_kutta(OrdDifEquation *ode, double x0, double y0, unsigned int M, DoubleVector &y, double h);
};

#endif // EXAMPLE2_H
