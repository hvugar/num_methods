#ifndef FIRST_ORDER_LINEAR_ODE_EX1_H
#define FIRST_ORDER_LINEAR_ODE_EX1_H

#include <ode/lode1o.h>
#include <utils/random.h>
#include <float.h>
#include <math.h>
#include "../border_global.h"

#define TIME_MAX 100
#define TIME_STEP 0.01

class FirstOrderLinearSample
{
public:
    double x(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const;
    double dt(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const;
    double d2t(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const;
    double d3t(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const;
    double d4t(const PointNodeODE &node, unsigned int r UNUSED_PARAM) const;
    //virtual double d5t(const PointNodeODE &node, unsigned int row = 0) const;
    //virtual double d6t(const PointNodeODE &node, unsigned int row = 0) const;

    void printNorms(std::vector<DoubleVector> &x, unsigned int k) const;
};

/*****************************************************************************************************/

class BORDERSHARED_EXPORT FirstOrderLinearODEEx1 : public IFirstOrderLinearODE, public FirstOrderLinearSample
{
public:
    FirstOrderLinearODEEx1();
    virtual ~FirstOrderLinearODEEx1();

    static void Main(int argc, char** argv);
    static void NonLocalConditionExample();
    static void CauchyProblemExample();

protected:
    virtual auto A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, unsigned int row = 0) const -> double;
    virtual auto count() const -> unsigned int;
    virtual auto dimension() const -> Dimension;

    virtual auto initial(InitialCondition, unsigned int) const -> double;
    virtual auto boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double { return 0.0; }
};

/*****************************************************************************************************/

class BORDERSHARED_EXPORT FirstOrderLinearODEFBVP : public IFirstOrderLinearODEFBVP, public FirstOrderLinearSample
{
protected:
    virtual auto A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, unsigned int row = 0) const -> double;
    virtual auto count() const -> unsigned int;
    virtual auto dimension() const -> Dimension;

    virtual auto final(FinalCondition, unsigned int) const -> double;
    virtual auto boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double { return 0.0; }
};

#endif // FIRST_ORDER_LINEAR_ODE_EX1_H
