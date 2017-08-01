#ifndef DIFFENSIALEQUATION_H
#define DIFFENSIALEQUATION_H

#include "global.h"
#include <vector2d.h>

class MINIMUMSHARED_EXPORT DiffensialEquation {};

class MINIMUMSHARED_EXPORT LinearODE : public DiffensialEquation {};

/**
 * @brief The Linear ODE 1st order in form dy/dx=a(x)y(x) + b(x);
 */
class MINIMUMSHARED_EXPORT LinearODE1stOrder : public LinearODE
{
protected:
    virtual double a(double t, unsigned int i) const = 0;
    virtual double b(double t, unsigned int i) const = 0;
};

/**
 * @brief The Linear ODE 2nd order in form d^2y/dx^2=q(x)dy/dx + p(x)*y(x) + r(x);
 */
class MINIMUMSHARED_EXPORT LinearODE2ndOrder : public LinearODE
{
protected:
    virtual double q(double x, unsigned int i) const = 0;
    virtual double p(double x, unsigned int i) const = 0;
    virtual double r(double x, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE : public DiffensialEquation {};

class MINIMUMSHARED_EXPORT NonLinearODE1stOrder : public NonLinearODE
{
protected:
    virtual double f(double x, double y, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE2ndOrder : public NonLinearODE {};

class MINIMUMSHARED_EXPORT SystenDiffensialEquation {};

class MINIMUMSHARED_EXPORT SystemLinearODE : public SystenDiffensialEquation {};

class MINIMUMSHARED_EXPORT SystemLinearODE1stOrder : public SystemLinearODE
{
protected:
    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t, unsigned int k, unsigned int row) const = 0;
};

class MINIMUMSHARED_EXPORT SystemNonLinearODE : public SystenDiffensialEquation {};

class MINIMUMSHARED_EXPORT SystemNonLinearODE1stOrder : public SystemNonLinearODE
{
protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const = 0;
};

#endif // DIFFENSIALEQUATION_H
