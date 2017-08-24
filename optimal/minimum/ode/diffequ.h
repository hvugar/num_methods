#ifndef DIFFENSIALEQUATION_H
#define DIFFENSIALEQUATION_H

#include "global.h"
#include <vector2d.h>
#include "grid/grid.h"

class MINIMUMSHARED_EXPORT DifferentialEquation
{
public:
    virtual unsigned int equationsNumber() const;

    const ODEGrid &grid() const;
    void setGrid(const ODEGrid& grid);
protected:
    ODEGrid mgrid;
};

class MINIMUMSHARED_EXPORT OrdinaryDifferentialEquation : public DifferentialEquation
{
public:
    enum Method
    {
        RK2,
        RK4,
        EULER,
        EULER_MOD
    };

    enum Direction
    {
        L2R, // Left to Right
        R2L  // Right to Left
    };
};

class MINIMUMSHARED_EXPORT ODE1stOrder {};

class MINIMUMSHARED_EXPORT ODE2ndOrder {};

class MINIMUMSHARED_EXPORT LinearODE : public OrdinaryDifferentialEquation {};

class MINIMUMSHARED_EXPORT NonLinearODE : public OrdinaryDifferentialEquation {};

/**
 * @brief The Linear ODE 2nd order in canonical (normal) form y"(x) = q(x)y'(x) + p(x)*y(x) + r(x);
 */
class MINIMUMSHARED_EXPORT LinearODE2ndOrder : public LinearODE, public ODE2ndOrder
{
protected:
    virtual double q(double x, unsigned int i) const = 0;
    virtual double p(double x, unsigned int i) const = 0;
    virtual double r(double x, unsigned int i) const = 0;
};

/**
 * @brief The NonLinear ODE1 2nd order in canonical (normal) form y"(x)=f(x,y(x),y'(x));
 */
class MINIMUMSHARED_EXPORT NonLinearODE2ndOrder : public NonLinearODE, public ODE2ndOrder
{
protected:
    virtual double f(double x, double y, double dy, unsigned int i) const = 0;
};

class MINIMUMSHARED_EXPORT SystemDifferentialEquation
{
public:
    const ODEGrid& grid() const;
    void setGrid(const ODEGrid& grid);
protected:
    ODEGrid mgrid;
};

class MINIMUMSHARED_EXPORT SystemDifferentialEquationODE : public SystemDifferentialEquation
{};

class MINIMUMSHARED_EXPORT SystemLinearODE : public SystemDifferentialEquationODE
{};

class MINIMUMSHARED_EXPORT SystemLinearODE1stOrder : public SystemLinearODE
{
protected:
    virtual double A(double t, unsigned int k, unsigned int row = 0, unsigned int col = 0) const = 0;
    virtual double B(double t, unsigned int k, unsigned int row = 0) const = 0;
};

class MINIMUMSHARED_EXPORT SystemNonLinearODE : public SystemDifferentialEquationODE
{};

class MINIMUMSHARED_EXPORT SystemNonLinearODE1stOrder : public SystemNonLinearODE
{
protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const = 0;
};

#endif // DIFFENSIALEQUATION_H
