#ifndef STEEPEST_DESCENT_GRADIENT_H
#define STEEPEST_DESCENT_GRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Steepest Descent Gradient
 */
class MINIMUMSHARED_EXPORT SteepestDescentGradient : public GradientMethod, protected R1Function
{
public:
    SteepestDescentGradient();
    virtual ~SteepestDescentGradient();

    virtual void calculate();
    virtual void calculate(DoubleVector &x);

protected:
    virtual double minimize();
    virtual double fx(double alpha);
};

#endif // STEEPEST_DESCENT_GRADIENT_H
