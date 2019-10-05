#ifndef FIRST_ORDER_NON_LINEAR_ODE_EX1_H
#define FIRST_ORDER_NON_LINEAR_ODE_EX1_H

#include <ode/nlode1o.h>
#include <utils/random.h>
#include <cfloat>
#include <cmath>

class MINIMUMSHARED_EXPORT FirstOrderNonLinearODErEx1 : public FirstOrderNonLinearODE
{
public:
    static void Main(int argc, char** argv);

protected:
    virtual unsigned int count() const;
    virtual double x(const PointNodeODE &node, unsigned int row = 1) const;
    virtual double f(const PointNodeODE &node, const DoubleVector &x, unsigned int r) const;
};

class FirstOrderNonLinearODEEx2 : public FirstOrderNonLinearODE
{
public:
    virtual double f(double x, double y, unsigned int k) const;
};

#endif // FIRST_ORDER_NON_LINEAR_ODE_EX1_H
