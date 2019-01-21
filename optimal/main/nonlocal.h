#ifndef NONLOCAL_H
#define NONLOCAL_H

#include <vector2d.h>
#include <limits>
#include <cmath>
#include <cfloat>
#include <linearequation.h>
#include <vector>

class NonLocal
{
public:
    NonLocal();
    virtual ~NonLocal();

    virtual void calculateAlphaSingle(double *alpha, unsigned int n, unsigned int k, double h) const;
    virtual void calculateAlphaSystem(DoubleMatrix *alpha, unsigned int n, unsigned int k, double h, unsigned int M) const;

    virtual void calculateOtherSingle(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, double *betta) const;
    virtual void calculateOtherSystem(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, DoubleMatrix *betta, unsigned int M) const;

    void solveSingle(const DoubleVector &a, double b, double h) const;
    void solveSystem(const std::vector<DoubleMatrix> &C, const DoubleVector &d, double h, unsigned int N) const;

    virtual double a(double t) const;
    virtual double b(double t) const;

    virtual double A(double t, unsigned int r, unsigned int c) const;
    virtual double B(double t, unsigned int n) const;

    virtual double x(double t) const;
    virtual double x(double t, unsigned int n) const;
};

#endif // NONLOCAL_H
