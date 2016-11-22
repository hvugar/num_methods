#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"
#include "printer.h"
#include "projection.h"

class ConjugateGradient : public GradientMethod, protected R1Function
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate(DoubleVector &x);
protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g);
    virtual double fx(double alpha);
    
    DoubleVector *mx;
    DoubleVector *ms;
};

#endif // CONJUGATEGRADIENT_H
