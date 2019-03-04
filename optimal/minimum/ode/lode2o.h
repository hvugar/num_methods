#ifndef SECOND_ORDER_LINEAR_ODE_H
#define SECOND_ORDER_LINEAR_ODE_H

#include "diffequ.h"
#include "cmethods.h"

class BoundaryConditionODE
{
public:
    BoundaryCondition boundaryConditionType;
    DoubleMatrix a;
    DoubleMatrix b;
    double lambda;
};

class InitialConditionODE
{
public:
    InitialCondition initialConditionType;
    double value;
};

/**
 * @brief Линейное дифференциальное уравнение второго порядка с переменными коэффициентами
 * The Linear ODE 2nd order in canonical (normal) form R(x)y"(x) + P(x)y'(x) + Q(x)y(x) = F(x);
 */
class MINIMUMSHARED_EXPORT SecondOrderLinearODE : public LinearODE
{
public:
    /**
     * @brief solveInitialValueProblem
     * @param rv
     */
    void solveInitialValueProblem(DoubleVector &rv) const;
    /**
     * @brief solveBoundaryValueProblem
     * @param rv
     */
    void solveBoundaryValueProblem(DoubleVector &rv) const;
    /**
     * @brief solveInitialValueProblem
     * @param rv
     */
    void solveInitialValueProblem(std::vector<DoubleVector> &rv) const;
    /**
     * @brief solveBoundaryValueProblem
     * @param rv
     */
    void solveBoundaryValueProblem(std::vector<DoubleVector> &rv) const;

protected:
    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;
    virtual auto B(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double = 0;

    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double = 0;
};

#endif // SECOND_ORDER_LINEAR_ODE_H
