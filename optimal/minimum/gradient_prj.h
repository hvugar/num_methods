#ifndef PRJGRADIENT_H
#define PRJGRADIENT_H

#include "gradient_sd.h"

class MINIMUMSHARED_EXPORT ProjectionGradient : public SteepestDescentGradient
{
public:
    ProjectionGradient();
    virtual ~ProjectionGradient();

    virtual void calculate();
    virtual void calculate(DoubleVector& x);

    double a;
    double b;

protected:
    double fx(double alpha);
};

#endif // PRJGRADIENT_H
