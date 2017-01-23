#ifndef LINEARBOUNDARYVALUEPROBLEMODE_H
#define LINEARBOUNDARYVALUEPROBLEMODE_H

#include "bvp.h"

/**
 * @brief The LinearBoundaryValueProblemODE class
 * r(x)y"(x)+p(x)y'(x)+q(x)y(x)=f(x);
 * y(0) = boundary(Left);
 * y(N) = boundart(Right);
 * @see BoundaryValueProblemODE
 * @see BoundaryValueProblem
 */
class MINIMUMSHARED_EXPORT LinearBoundaryValueProblemODE : protected BoundaryValueProblemODE
{
protected:
    virtual double r(unsigned int n) const = 0;
    virtual double p(unsigned int n) const = 0;
    virtual double q(unsigned int n) const = 0;
    virtual double f(unsigned int n) const = 0;

public:
    void calculate(DoubleVector &x, unsigned int k, double h, unsigned int N, unsigned int direction = 1);
    void calculate2N(DoubleVector &x, double h, unsigned int N);

    void calculateN4L2RD(DoubleVector &x, double h, unsigned int N);
    void calculateN4R2LD(DoubleVector &x, double h, unsigned int N);
    void calculateN6L2RD(DoubleVector &x, double h, unsigned int N);
    void calculateN6R2LD(DoubleVector &x, double h, unsigned int N);
};

#endif // LINEARBOUNDARYVALUEPROBLEMODE_H
