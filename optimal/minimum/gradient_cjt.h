#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Conjugate Gradient
 */
class MINIMUMSHARED_EXPORT ConjugateGradient : public GradientMethod
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate(DoubleVector& x);
protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g);
};

#endif // CONJUGATEGRADIENT_H
