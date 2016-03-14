#ifndef BORDERPARABOLIC1D321_H
#define BORDERPARABOLIC1D321_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <gradient.h>
#include <function.h>
#include <math.h>

class Parabolic1DControl321 : public IParabolicEquation, public RnFunction
{
public:
    Parabolic1DControl321();
    virtual ~Parabolic1DControl321() {}

    virtual double fx(const DoubleVector &x);

    virtual double fi(unsigned int i) const;
    virtual double m1(unsigned int j) const;
    virtual double m2(unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

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

    DoubleVector e;

    double v1(double t) const;
    double v2(double t) const;

    DoubleVector U;
};

#endif // BORDERPARABOLIC1D321_H
