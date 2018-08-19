#ifndef STEEPEST_DESCENT_GRADIENT_H
#define STEEPEST_DESCENT_GRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Steepest Descent Gradient
 * Метод наискорейшего градиентного спуска.
 */
class MINIMUMSHARED_EXPORT SteepestDescentGradient : public GradientMethod, protected R1Function
{
public:
    SteepestDescentGradient();
    virtual ~SteepestDescentGradient();

    virtual void calculate(DoubleVector &x);

    R1FxMinimizer &r1Minimizer();
    const R1FxMinimizer& r1Minimizer() const;

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;
    virtual double fx(double alpha) const;

    DoubleVector *mx;
    DoubleVector *mg;

    R1FxMinimizer r1m;
};

#endif // STEEPEST_DESCENT_GRADIENT_H
