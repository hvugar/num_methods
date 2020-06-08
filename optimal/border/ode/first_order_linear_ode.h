#ifndef FIRST_ORDER_LINEAR_ODE_EX1_H
#define FIRST_ORDER_LINEAR_ODE_EX1_H

#include "../border_global.h"

#define TIME_MAX 100
#define TIME_STEP 0.01

class FirstOrderLinearSample
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

class BORDERSHARED_EXPORT FirstOrderLinearODEIVP : public IFirstOrderLinearODEIVP, public FirstOrderLinearSample
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();
    static void NonLocalConditionExample();

    std::vector<DoubleVector> mx;

protected:
    virtual auto A(const PointNodeODE &node, size_t row = 1, size_t col = 1) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 1) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual void iterationInfo(double y, const PointNodeODE &node) const;
    virtual void iterationInfo(const DoubleVector &v, const PointNodeODE &node) const;
    virtual auto initial(InitialCondition, size_t) const -> double;
};

/*****************************************************************************************************/

class BORDERSHARED_EXPORT FirstOrderLinearODEFVP : public IFirstOrderLinearODEFVP, public FirstOrderLinearSample
{
public:
    static void Main(int argc, char** argv);
    static void CauchyProblemExample();

    std::vector<DoubleVector> mx;

protected:
    virtual auto A(const PointNodeODE &node, size_t row = 0, size_t col = 0) const -> double;
    virtual auto B(const PointNodeODE &node, size_t row = 0) const -> double;
    virtual auto count() const -> size_t;
    virtual auto dimension() const -> Dimension;
    virtual auto iterationInfo(double y, const PointNodeODE &node) const -> void;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void;
    virtual auto final(FinalCondition, size_t) const -> double;
};

#endif // FIRST_ORDER_LINEAR_ODE_EX1_H
