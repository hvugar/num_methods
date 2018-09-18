#ifndef STEEPEST_DESCENT_GRADIENT_H
#define STEEPEST_DESCENT_GRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Steepest Descent Gradient
 * Метод наискорейшего градиентного спуска.
 */
class MINIMUMSHARED_EXPORT SteepestDescentGradient : public GradientMethod
{
public:
    SteepestDescentGradient();
    virtual ~SteepestDescentGradient();

    virtual void calculate(DoubleVector &x);

//protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;
    virtual double fx(double alpha) const;

    DoubleVector *mx;
    DoubleVector *mg;
};

#endif // STEEPEST_DESCENT_GRADIENT_H
