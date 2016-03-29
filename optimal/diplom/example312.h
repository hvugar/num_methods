#ifndef BORDERPARABOLIC1D312_H
#define BORDERPARABOLIC1D312_H

#include <parabolicequation.h>
#include <doublevector.h>
#include <math.h>

class BorderParabolic1D312 : public IParabolicEquation
{
public:
    BorderParabolic1D312();
    virtual ~BorderParabolic1D312() {}

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
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

    DoubleVector e;
    double e1;
    double e2;
    double e3;

    double v1(double t) const;
    double v2(double t) const;
    double v3(double t) const;
};

#endif // BORDERPARABOLIC1D2_H
