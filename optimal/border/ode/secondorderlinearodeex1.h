#ifndef SECOND_ORDER_LINEAR_ODE_EX1_H
#define SECOND_ORDER_LINEAR_ODE_EX1_H

#include <ode/lode2o.h>

class MINIMUMSHARED_EXPORT SecondOrderLinearODEEx1 : public SecondOrderLinearODE
{
public:
    static void Main(int argc, char **argv);

protected:
    virtual double A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const;
    virtual double B(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const;
    virtual double C(const PointNodeODE &node, unsigned int row = 1) const;

    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double;
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double;

    virtual unsigned int count() const;
};

#endif // SECOND_ORDER_LINEAR_ODE_EX1_H
