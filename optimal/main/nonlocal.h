#ifndef NONLOCAL_H
#define NONLOCAL_H

#include <vector2d.h>
#include <limits>
#include <cmath>
#include <cfloat>
#include <linearequation.h>
#include <vector>
#include <grid/grid.h>
#include <ode/diffequ.h>
#include <utils/random.h>

class INonLocal
{
public:
    virtual ~INonLocal();

    virtual void calculateAlphaSingle(double *alpha, unsigned int n, unsigned int k, double h) const;
    virtual void calculateOtherSingle(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, double *betta) const;
    virtual void solveSingle(const DoubleVector &a, double b, double h) const;

protected:
    virtual double a(double t) const = 0;
    virtual double b(double t) const  = 0;
    virtual double x(double t) const = 0;
};

class INonLocalSystem : public SystemLinearODE1stOrder
{
public:
    virtual ~INonLocalSystem();

    virtual void solve(const std::vector<DoubleMatrix> &C, const DoubleVector &d, std::vector<DoubleVector> &x, const Dimension &dim, unsigned int N) const;

protected:
    virtual double A(double t, unsigned int n, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t, unsigned int n, unsigned int m) const = 0;

private:
    virtual void getAlpha(DoubleMatrix *alpha, unsigned int n, unsigned int k, double h, unsigned int M) const;
    virtual void getEquations(DoubleMatrix &c, DoubleVector &d, unsigned int N, unsigned int k, double h, DoubleMatrix *betta, unsigned int M) const;
};

////////////// Concrete classes ///////////////////////////////////////

class NonLocal : public INonLocal
{
public:
    virtual ~NonLocal();

    virtual double a(double t) const;
    virtual double b(double t) const;
    virtual double x(double t) const;
};

class NonLocalSystem : public INonLocalSystem
{
public:
    static void Main(int argc, char** argv);
    virtual ~NonLocalSystem();

    virtual double A(double t, unsigned int n, unsigned int row, unsigned int col) const;
    virtual double B(double t, unsigned int n, unsigned int m) const;
    virtual double x(double t, unsigned int n) const;
};

#endif // NONLOCAL_H
