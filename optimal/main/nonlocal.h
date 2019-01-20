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

    void solve(const DoubleVector &a, double b, double h) const;
    virtual void calculateAlpha(double *alpha, unsigned int n, unsigned int k, double h) const;
    virtual void calculateAlpha(DoubleVector *alpha, unsigned int m, unsigned int k, double h, unsigned int N) const;
    virtual void calculateOther(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, double *betta) const;

    void solveSystem(const std::vector<DoubleMatrix> &C, const DoubleVector &d, double h, unsigned int N) const;

    virtual double a(double t) const;
    virtual double b(double t) const;

    virtual double A(double t, unsigned int n) const;
    virtual double B(double t, unsigned int n) const;

    bool equalZero(double) const;

    virtual double x(double t) const;
};

#endif // NONLOCAL_H
