#ifndef CONSTSTEPGRADIENT_H
#define CONSTSTEPGRADIENT_H

#include "gradient.h"

/**
 * @brief Метод градиентного спуска с постоянным шагом.
 *
 */
class MINIMUMSHARED_EXPORT ConstStepGradient : public GradientBasedMethod
{
public:
    ConstStepGradient();
    virtual ~ConstStepGradient();

    virtual void calculate(DoubleVector &x);

protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g) const;
    virtual double fx(double alpha) const;

private:
    DoubleVector *mx;
    DoubleVector *mg;
};

#endif // CONSTSTEPGRADIENT_H
