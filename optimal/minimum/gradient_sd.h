#ifndef STEEPEST_DESCENT_GRADIENT_H
#define STEEPEST_DESCENT_GRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Steepest Descent Gradient
 */
class MINIMUMSHARED_EXPORT SteepestDescentGradient : public Gradient, protected R1Function
{
public:
    SteepestDescentGradient();
    virtual ~SteepestDescentGradient();

    virtual void calculate();

protected:
    virtual double minimize();
    virtual double fx(double alpha);
    virtual void print();
};

#endif // STEEPEST_DESCENT_GRADIENT_H
