#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "gradient.h"
#include "printer.h"
#include "projection.h"

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
    virtual double fx(double alpha);
    
    DoubleVector *mx;
    DoubleVector *ms;
};

class MINIMUMSHARED_EXPORT IConjugateGradient : public GradientMethod, public IGradient, public RnFunction, protected R1Function, public IProjection, public IPrinter
{
public:
    IConjugateGradient() {}
    virtual ~IConjugateGradient() = 0;

    virtual void gradient(const DoubleVector &x, DoubleVector &g) = 0;
    virtual double fx(const DoubleVector &x) = 0;
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &gradient, double alpha, RnFunction *fn) const = 0;

    virtual void calculate(DoubleVector &x);
protected:
    virtual double fx(double x) = 0;
    virtual double minimize(const DoubleVector &x, const DoubleVector &g);

private:
    DoubleVector *mx;
    DoubleVector *ms;
};

#endif // CONJUGATEGRADIENT_H
