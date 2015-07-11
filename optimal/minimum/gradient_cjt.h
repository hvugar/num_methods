#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

class MINIMUMSHARED_EXPORT ConjugateGradient : public Gradient, protected R1Function
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void setX(const std::vector<double>& x);

    virtual void calculate();

protected:
    virtual double minimize();
    virtual double distance() const;
    virtual double fx(double alpha);
    virtual void print();

private:
    std::vector<double> s;
};

#endif // CONJUGATEGRADIENT_H
