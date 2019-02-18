#ifndef SECOND_ORDER_NON_LINEAR_ODE_H
#define SECOND_ORDER_NON_LINEAR_ODE_H

#include "diffequ.h"

/**
 * @brief The NonLinear ODE1 2nd order in canonical (normal) form y"(x)=f(x,y(x),y'(x));
 */
class SecondOrderNonLinearODE : virtual public NonLinearODE
{
protected:
    virtual double f(double x, double y, double dy, unsigned int k) const;
    virtual double f(double x, const DoubleVector &y, const DoubleVector &dy, unsigned int k, unsigned int i) const;
};

#endif // SECOND_ORDER_NON_LINEAR_ODE_H
