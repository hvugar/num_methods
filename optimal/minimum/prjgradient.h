#ifndef PRJGRADIENT_H
#define PRJGRADIENT_H

#include "gradient.h"
#include "sdgradient.h"

class MINIMUMSHARED_EXPORT ProjectionGradient : public SteepestDescentGradient
{
public:
    ProjectionGradient();
    virtual ~ProjectionGradient();

    virtual void calculate();

    double a;
    double b;

protected:
    double fx(double alpha);
};

#endif // PRJGRADIENT_H
