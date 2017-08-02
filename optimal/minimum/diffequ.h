#ifndef DIFFENSIALEQUATION_H
#define DIFFENSIALEQUATION_H

#include "global.h"
#include <vector2d.h>
#include "grid/grid.h"

class MINIMUMSHARED_EXPORT DifferentialEquation {};

class MINIMUMSHARED_EXPORT OrdinaryDifferentialEquation : protected DifferentialEquation {};

class MINIMUMSHARED_EXPORT ODE1stOrder {};

class MINIMUMSHARED_EXPORT ODE2ndOrder {};

class MINIMUMSHARED_EXPORT LinearODE : public OrdinaryDifferentialEquation {};

/**
 * @brief The Linear ODE 1st order in form dy/dx=a(x)y(x) + b(x);
 */
class MINIMUMSHARED_EXPORT LinearODE1stOrder : public LinearODE, public ODE1stOrder
{
protected:
    virtual double a(double t, unsigned int i) const = 0;
    virtual double b(double t, unsigned int i) const = 0;
};

/**
 * @brief The Linear ODE 2nd order in form d^2y/dx^2=q(x)dy/dx + p(x)*y(x) + r(x);
 */
class MINIMUMSHARED_EXPORT LinearODE2ndOrder : public LinearODE, public ODE2ndOrder
{
protected:
    virtual double q(double x, unsigned int i) const = 0;
    virtual double p(double x, unsigned int i) const = 0;
    virtual double r(double x, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE : public OrdinaryDifferentialEquation {};

class MINIMUMSHARED_EXPORT NonLinearODE1stOrder : public NonLinearODE, public ODE1stOrder
{
protected:
    virtual double f(double x, double y, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT NonLinearODE2ndOrder : public NonLinearODE, public ODE2ndOrder {};



class MINIMUMSHARED_EXPORT SystemDifferentialEquation
{
public:
    SystemDifferentialEquation(const ODEGrid &grid);
    const ODEGrid& grid() const;
protected:
    ODEGrid mgrid;
};

class MINIMUMSHARED_EXPORT SystemDifferentialEquationODE : public SystemDifferentialEquation
{
public:
    SystemDifferentialEquationODE(const ODEGrid &grid) : SystemDifferentialEquation(grid) {}
};

class MINIMUMSHARED_EXPORT SystemLinearODE : public SystemDifferentialEquationODE
{
public:
    SystemLinearODE(const ODEGrid &grid) : SystemDifferentialEquationODE(grid) {}
};

class MINIMUMSHARED_EXPORT SystemLinearODE1stOrder : public SystemLinearODE
{
public:
    SystemLinearODE1stOrder(const ODEGrid &grid);
protected:
    virtual double A(double t, unsigned int k, unsigned int row, unsigned int col) const = 0;
    virtual double B(double t, unsigned int k, unsigned int row) const = 0;
};

class MINIMUMSHARED_EXPORT SystemNonLinearODE : public SystemDifferentialEquationODE
{
public:
    SystemNonLinearODE(const ODEGrid &grid) : SystemDifferentialEquationODE(grid) {}
};

class MINIMUMSHARED_EXPORT SystemNonLinearODE1stOrder : public SystemNonLinearODE
{
public:
    SystemNonLinearODE1stOrder(const ODEGrid &grid) : SystemNonLinearODE(grid) {}
protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const = 0;
};

#endif // DIFFENSIALEQUATION_H
