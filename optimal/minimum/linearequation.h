#ifndef LINEAR_EQUATION_H
#define LINEAR_EQUATION_H

#include "global.h"
#include "matrix2d.h"

//void MINIMUMSHARED_EXPORT GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector &x);

class MINIMUMSHARED_EXPORT LinearEquation
{
public:
    static void GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector& x);
};

#endif // LINEAR_EQUATION_H
