#ifndef BORDERPARABOLIC1D1_H
#define BORDERPARABOLIC1D1_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <math.h>

class BorderParabolic1D1 : public IParabolicEquation
{
public:
    BorderParabolic1D1();
    virtual ~BorderParabolic1D1() {}

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

    double e1;
    double e2;
    double e3;

    double v1(double t) const;
    double v2(double t) const;
};

#endif // BORDERPARABOLIC1D1_H
