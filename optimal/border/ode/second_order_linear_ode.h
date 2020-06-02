#ifndef SECOND_ORDER_LINEAR_ODE_EX1_H
#define SECOND_ORDER_LINEAR_ODE_EX1_H

#include "../border_global.h"
#include <ode/lode2o.h>

class BORDERSHARED_EXPORT SecondOrderLinearODEIBVP : public ISecondOrderLinearODEIBVP
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();

    virtual auto A(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto C(const PointNodeODE &node, size_t row = 0) const -> double;
    virtual auto initial(InitialCondition, size_t) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;

    virtual auto boundary(const PointNodeODE &node, BoundaryConditionPDE &condition, size_t row = 1) const -> double;
};

class BORDERSHARED_EXPORT SecondOrderLinearODEFBVP : public ISecondOrderLinearODEFBVP
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();

    virtual auto A(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto C(const PointNodeODE &node, size_t row = 0) const -> double;
    virtual auto final(FinalCondition, size_t) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;

    virtual auto boundary(const PointNodeODE &node, BoundaryConditionPDE &condition, size_t row = 1) const -> double;
};

#endif // SECOND_ORDER_LINEAR_ODE_EX1_H
