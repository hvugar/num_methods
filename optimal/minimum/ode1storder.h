#ifndef ODE1STORDER_H
#define ODE1STORDER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cmethods.h"
#include "global.h"

class MINIMUMSHARED_EXPORT DifferensialEquation {};

class MINIMUMSHARED_EXPORT LinearODE : public DifferensialEquation {};

/**
 * @brief The Linear ODE 1st order in form dy/dx=a(x)y(x) + b(x);
 */
class MINIMUMSHARED_EXPORT LinearODE_1stOrder : public LinearODE
{
protected:
    virtual double a(double t, unsigned int i) const = 0;
    virtual double b(double t, unsigned int i) const = 0;
};

/**
 * @brief The Linear ODE 2nd order in form d^2y/dx^2=q(x)dy/dx + p(x)*y(x) + r(x);
 */
class MINIMUMSHARED_EXPORT LinearODE_2ndOrder : public LinearODE
{
protected:
    virtual double q(double x, unsigned int i) const = 0;
    virtual double p(double x, unsigned int i) const = 0;
    virtual double r(double x, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE : public DifferensialEquation {};

class MINIMUMSHARED_EXPORT NonLinearODE1Order : public NonLinearODE
{
protected:
    virtual double f(double x, double y, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE2ndOrder : public NonLinearODE {};

class MINIMUMSHARED_EXPORT LinearODE1stOrder
{
public:
    virtual double a(double t, unsigned int i) const = 0;
    virtual double b(double t, unsigned int i) const = 0;
    virtual double initial(double t, unsigned int i) const = 0;
    virtual double boundary(double t, unsigned int i) const = 0;
public:
    void solveLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const;
};

class MINIMUMSHARED_EXPORT NonLinearODE1stOrder
{
    virtual double f(double t, double x, unsigned int i) const = 0;
    virtual double initial(double t, unsigned int i) const = 0;
    virtual double boundary(double t, unsigned int i) const = 0;
public:
    void solveNonLinearBoundaryProblem(double t0, double x0, double tn, double xn, unsigned int N, double *x) const;
};

#endif // ODE1STORDER_H
