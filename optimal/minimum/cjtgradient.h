#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT ConjugateGradient : public Gradient
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    void calculate();
    double minimize();
    void print();

private:
    std::vector<double> s;
};

#endif // CONJUGATEGRADIENT_H
