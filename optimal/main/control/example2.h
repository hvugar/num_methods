#ifndef EXAMPLE2_H
#define EXAMPLE2_H

#include <function.h>
#include <printer.h>
#include <stdlib.h>
#include <float.h>
#include <function.h>
#include <gradient.h>
#include <gradient_cjt.h>

class Ex2OrdDifEquationX;
class Ex2OrdDifEquationP;
class Example2;

struct Ex2OrdDifEquationX
{
    Example2 *e;
    virtual double fx(double t, double x) const;
};

struct Ex2OrdDifEquationP
{
    Example2 *e;
    virtual double fx(double t, double p, unsigned int i) const;
};

class Example2 : public ConjugateGradient, public RnFunction, public IGradient, public IPrinter
{
public:
    Example2();

    virtual double fx(const DoubleVector &x);
    virtual void gradient(const DoubleVector &x, DoubleVector &g);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const;

    static void main(int argc, char* argv[]);

    void calculateX(DoubleVector &x);
    void calculateP(DoubleVector &p);

    double A(double t) { return 3.0*t; }
    double B(double t) { return 1.0; }
    double C(double t) { return 2.0*t - 3.0*t*t; }

private:
    double t0;
    double t1;
    double h;
    unsigned int M;
    // Number of t interval seperators; Interval count is N+1. N+1 ~ T
    unsigned int N;
    // Number of x interval seperators; Interval count is S+1
    unsigned int S;
    DoubleVector xs;


    double x0;
    double xT;
    double xi;

    const DoubleVector *pK;
    const DoubleVector *px;
    const DoubleVector *pp;

    double getK(double x) const;
    unsigned int getKS(double x) const;

    Ex2OrdDifEquationX eqX;
    Ex2OrdDifEquationP eqP;

    friend class Ex2OrdDifEquationX;
    friend class Ex2OrdDifEquationP;
};

#endif // EXAMPLE2_H
