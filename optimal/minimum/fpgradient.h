#ifndef FPGRADIENT_H
#define FPGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT FastProximalGradient : public Gradient
{
public:
    FastProximalGradient();
    virtual ~FastProximalGradient();

    void calculate();
    void print();
};

#endif // FPGRADIENT_H
