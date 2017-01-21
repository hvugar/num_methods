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
    void calculate2N(DoubleVector &x, double h, unsigned int N);
    void calculate4NL2R(DoubleVector &x, double h, unsigned int N);
    void calculate4NR2L(DoubleVector &x, double h, unsigned int N);
    void calculate6NL2R(DoubleVector &x, double h, unsigned int N);
    void calculate6NR2L(DoubleVector &x, double h, unsigned int N);
};

#endif // LINEARBOUNDARYVALUEPROBLEMODE_H
