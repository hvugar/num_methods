#ifndef BOUNDARYVALUEPROBLEM_H
#define BOUNDARYVALUEPROBLEM_H

#include <vector2d.h>
#include <matrix2d.h>
#include <printer.h>

/**
 * @brief The BoundaryProblem class
 * r(x)y"(x)+p(x)y'(x)+q(x)y(x)=f(x);
 * y(0) = boundary(Left);
 * y(N) = boundart(Right);
 */

class MINIMUMSHARED_EXPORT BoundaryValueProblem
{
protected:
    virtual double r(unsigned int i) const = 0;
    virtual double p(unsigned int i) const = 0;
    virtual double q(unsigned int i) const = 0;
    virtual double f(unsigned int i) const = 0;
    virtual double boundary(Boundary bound) const = 0;
public:
    void calculate2N(DoubleVector &x, double h, unsigned int N);
    void calculate4NL2R(DoubleVector &x, double h, unsigned int N);
    void calculate4NR2L(DoubleVector &x, double h, unsigned int N);
    void calculate6NL2R(DoubleVector &x, double h, unsigned int N);
    void calculate6NR2L(DoubleVector &x, double h, unsigned int N);
};

#endif // BOUNDARYVALUEPROBLEM_H
