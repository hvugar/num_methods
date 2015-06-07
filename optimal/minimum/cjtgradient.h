#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT ConjugateGradient : public Gradient, protected R1Function
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate();
    virtual double minimize();
    virtual double distance() const;
    virtual void print();

protected:
    virtual double fx(double alpha);
    std::vector<double> s;
};

#endif // CONJUGATEGRADIENT_H
