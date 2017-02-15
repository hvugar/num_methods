#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"

/**
 * @brief Method of Conjugate Gradient
 * Метод Флетчера-Ривса (Метод сопряженных градиентов).
 */
class MINIMUMSHARED_EXPORT ConjugateGradient : public GradientMethod, protected R1Function
{
public:
    ConjugateGradient();
    virtual ~ConjugateGradient();

    virtual void calculate(DoubleVector &x);
protected:
    virtual double minimize(const DoubleVector &x, const DoubleVector &g);
    virtual double fx(double alpha) const;
    
    DoubleVector *mx;
    DoubleVector *ms;
};

#endif // CONJUGATEGRADIENT_H
