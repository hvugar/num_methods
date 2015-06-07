#ifndef FPGRADIENT_H
#define FPGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT FastProximalGradient : public Gradient, protected R1Function
{
public:
    FastProximalGradient();
    virtual ~FastProximalGradient();

    virtual double minimize();
    virtual void calculate();
    virtual void print();

protected:
    virtual double fx(double alpha);
};

#endif // FPGRADIENT_H
