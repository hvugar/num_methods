#ifndef LINEAR_EQUATION_H
#define LINEAR_EQUATION_H

#include "global.h"
#include "matrix2d.h"

//void MINIMUMSHARED_EXPORT GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector &x);

class MINIMUMSHARED_EXPORT LinearEquation
{
public:
    static void GaussianElimination(const DoubleMatrix& m, const DoubleVector& b, DoubleVector& x);
    static void FirstRowLoaded(const double *e, double f, const double *a, const double *b, const double *c, const double *d, unsigned int N);
};

#endif // LINEAR_EQUATION_H
