#ifndef SECOND_ORDER_LINEAR_ODE_H
#define SECOND_ORDER_LINEAR_ODE_H

#include "diffequ.h"
#include "cmethods.h"

enum class BoundaryConditionType
{
    Dirichlet = 1,
    Neumann = 2,
    Robin = 3
};

enum class InitialConditionType
{
    InitialValue = 0,
    InitialDerivative = 1
};

class BoundaryConditionODE
{
public:
    BoundaryConditionType boundaryConditionType;
    DoubleMatrix a;
    DoubleMatrix b;
    DoubleVector c;
};

class InitialConditionODE
{
public:
    InitialConditionType initialConditionType;
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
     * @brief solveBoundaryValueProblem
     * @param left  Boundary condition on left side
     * @param right Boundary condition on right side
     * @param ry    Result functions
     */
    void solveBoundaryValueProblem(const DoubleVector &left, const DoubleVector &right, std::vector<DoubleVector> &ry) const;
    /**
     * @brief solveBoundaryValueProblem
     * @param rv
     */
    void solveBoundaryValueProblem(std::vector<DoubleVector> &rv) const;
    /**
     * @brief solveBoundaryValueProblem
     * @param rv
     */
    void solveBoundaryValueProblem(DoubleVector &rv) const;
    /**
     * @brief solveInitialValueProblem
     * @param rv
     */
    void solveInitialValueProblem(DoubleVector &rv) const;
    /**
     * @brief solveCauchyProblem
     * @param initial Initial conditions
     * @param initder Initial derivative conditions
     * @param ry      Result functions
     */
    void solveCauchyProblem(const DoubleVector &initial, const DoubleVector &initder, std::vector<DoubleVector> &ry) const;

protected:
    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;
    virtual auto B(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double = 0;

    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition) const -> void = 0;
    virtual auto initial(const PointNodeODE &node, InitialConditionODE &condition) const -> void = 0;
};

#endif // SECOND_ORDER_LINEAR_ODE_H
