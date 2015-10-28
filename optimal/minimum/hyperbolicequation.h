#ifndef HYPERBOLICEQUATION_H
#define HYPERBOLICEQUATION_H

#include "global.h"
#include "doublevector.h"

class MINIMUMSHARED_EXPORT HyperbolicEquation
{
public:
    HyperbolicEquation(double t0 = 0.0, double t1 = 1.0, double x0 = 0.0, double x1 = 1.0, double a = 1.0, unsigned int M = 100, unsigned int N = 100);
    virtual ~HyperbolicEquation();

    void setTimeInterval(double t0, double t1);
    void setLengthInterval(double x0, double x1);
    void setPartNumber(unsigned int M, unsigned int N);

    virtual double fi(double x) const = 0;
    virtual double m1(double t) const = 0;
    virtual double m2(double t) const = 0;
    virtual double f(double x, double t) const = 0;

    virtual void calculate(DoubleVector& u);

private:
    double t0;
    double t1;
    double x0;
    double x1;
    double a;
    /* number of parts of time */
    unsigned int M;
    /* number of parts of length */
    unsigned int N;
    double ht;
    double hx;
};

#endif // HYPERBOLICEQUATION_H
