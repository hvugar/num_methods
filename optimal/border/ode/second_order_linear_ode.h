#ifndef SECOND_ORDER_LINEAR_ODE_EX1_H
#define SECOND_ORDER_LINEAR_ODE_EX1_H

#include "../border_global.h"

#define TIME_MIN_2ND  0
#define TIME_MAX_2ND  100
#define TIME_STEP_2ND 0.01

class SecondOrderLinearSample
{
public:
    double x(const PointNodeODE &node, size_t r UNUSED_PARAM) const;
    double dt(const PointNodeODE &node, size_t r UNUSED_PARAM) const;
    double d2t(const PointNodeODE &node, size_t r UNUSED_PARAM) const;
    double d3t(const PointNodeODE &node, size_t r UNUSED_PARAM) const;
    double d4t(const PointNodeODE &node, size_t r UNUSED_PARAM) const;
    //virtual double d5t(const PointNodeODE &node, size_t row = 0) const;
    //virtual double d6t(const PointNodeODE &node, size_t row = 0) const;

    void printNorms(std::vector<DoubleVector> &x, size_t k) const;
};

/*****************************************************************************************************/

class BORDERSHARED_EXPORT SecondOrderLinearODEIBVP : public ISecondOrderLinearODEIVP, public SecondOrderLinearSample
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();

    std::vector<DoubleVector> mx;

protected:
    virtual auto A(const PointNodeODE &node, size_t row = 1, size_t col = 1) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 1, size_t col = 1) const -> double;
    virtual auto C(const PointNodeODE &node, size_t row = 1) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;
    virtual auto initial(InitialCondition, size_t) const -> double;

    //virtual auto boundary(const PointNodeODE &node, BoundaryConditionPDE &condition, size_t row = 1) const -> double;
};

/*****************************************************************************************************/

class BORDERSHARED_EXPORT SecondOrderLinearODEFBVP : public ISecondOrderLinearODEFVP, public SecondOrderLinearSample
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();

    std::vector<DoubleVector> mx;

protected:
    virtual auto A(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto C(const PointNodeODE &node, size_t row = 0) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;
    virtual auto final(FinalCondition, size_t) const -> double;

    //virtual auto boundary(const PointNodeODE &node, BoundaryConditionPDE &condition, size_t row = 1) const -> double;
};

#endif // SECOND_ORDER_LINEAR_ODE_EX1_H
