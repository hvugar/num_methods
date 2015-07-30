#ifndef ROSENBROCK_H
#define ROSENBROCK_H

#include "function.h"

class Rosenbrock : public RnFunction
{
public:
    virtual double fx(const DoubleVector& x) const;
    virtual void gradient(DoubleVector& g) const;

    static void Main();

private:
    double grad_step;
};

#endif // ROSENBROCK_H
