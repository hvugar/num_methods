#ifndef FIRST_ORDER_LINEAR_ODE_EX1_H
#define FIRST_ORDER_LINEAR_ODE_EX1_H

#include <ode/lode1o.h>
#include <utils/random.h>
#include <float.h>
#include <math.h>
#include "../border_global.h"

class BORDERSHARED_EXPORT FirstOrderLinearODEEx1 : public FirstOrderLinearODE
{
public:
    FirstOrderLinearODEEx1();
    virtual ~FirstOrderLinearODEEx1();

    static void Main(int argc, char** argv);

protected:
    virtual double A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const;
    virtual double B(const PointNodeODE &node, unsigned int row = 0) const;
    virtual unsigned int count() const;

    virtual auto initial(InitialCondition, unsigned int) const -> double { return 0.0; }
    virtual auto boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double { return 0.0; }

    virtual double x(const PointNodeODE &node, unsigned int row = 0) const;
};

#endif // FIRST_ORDER_LINEAR_ODE_EX1_H
