#ifndef ODE1STORDER_H
#define ODE1STORDER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cmethods.h"
#include "global.h"

class MINIMUMSHARED_EXPORT LinearODE1stOrderIBVP
{
public:
    virtual double a(double t, unsigned int i) const = 0;
    virtual double b(double t, unsigned int i) const = 0;
    virtual double initial(double t, unsigned int i) const = 0;
    virtual double boundary(double t, unsigned int i) const = 0;
public:
    void solveLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const;
};

class MINIMUMSHARED_EXPORT NonLinearODE1stOrderIBVP
{
    virtual double f(double t, double x, unsigned int i) const = 0;
    virtual double initial(double t, unsigned int i) const = 0;
    virtual double boundary(double t, unsigned int i) const = 0;
public:
    void solveNonLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const;
};

#endif // ODE1STORDER_H
