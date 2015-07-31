#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Conjugate Gradient
 */
class MINIMUMSHARED_EXPORT ConjugateGradient : public GradientMethod, protected R1Function
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate();

    virtual void setX(const DoubleVector& x);

protected:
    virtual double minimize();
    virtual double distance() const;
    virtual double fx(double alpha);

private:
    DoubleVector s;
};

#endif // CONJUGATEGRADIENT_H
