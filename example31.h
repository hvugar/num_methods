#ifndef BORDERPARABOLIC1D311_H
#define BORDERPARABOLIC1D311_H

#include "parabolicequation.h"
#include "printer.h"
#include <math.h>

class BorderParabolic1D31 : public IParabolicEquation
{
public:
    BorderParabolic1D31();
    virtual ~BorderParabolic1D31() {}

    virtual double initial(unsigned int i) const;
    virtual double boundary(Boundary type, unsigned int j) const;
    virtual double f(unsigned int i, unsigned int j) const;

    static void main(int argc, char* argv[]);

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

#endif // BORDERPARABOLIC1D311_H
