#ifndef HYPERBOLICIBVP_H
#define HYPERBOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The HyperbolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t)
 * u(x,0) = fi1(x),
 * u_t(x,0) = fi2(x),
 * u(0,t) = m1(t),
 * u(1,t) = m2(t).
 */
class HyperbolicIBVP : protected InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial1(unsigned int n) const = 0;
    virtual double initial2(unsigned int n) const = 0;
    virtual double boundary(unsigned int m, BoundaryType boundary) const = 0;
    virtual double f(unsigned int n, unsigned int m) const = 0;
    virtual double a(unsigned int n, unsigned int m) const = 0;

public:
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
};

#endif // HYPERBOLICIBVP_H
