#ifndef NONLINEARODE2NDORDER_H
#define NONLINEARODE2NDORDER_H

#include "diffequ.h"

/**
 * @brief The NonLinear ODE1 2nd order in canonical (normal) form y"(x)=f(x,y(x),y'(x));
 */
class NonLinearODE2ndOrder : virtual public NonLinearODE
{
protected:
    virtual double f(double x, double y, double dy, unsigned int k) const;
    virtual double f(double x, const DoubleVector &y, const DoubleVector &dy, unsigned int k, unsigned int i) const;
};

#endif // NONLINEARODE2NDORDER_H
