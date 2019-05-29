#ifndef SECOND_ORDER_LINEAR_ODE_H
#define SECOND_ORDER_LINEAR_ODE_H

#include "diffequ.h"
#include "../grid/ibvp.h"

/**
 * @brief Линейное дифференциальное уравнение второго порядка с переменными коэффициентами
 * The Linear ODE 2nd order in canonical (normal) form y"(x) = A(x)y'(x) + B(x)y(x) + C(x);
 */
class MINIMUMSHARED_EXPORT SecondOrderLinearODE :
        virtual public LinearODE,
        virtual public InitialValueProblemODE,
        virtual public BoundaryValueProblemODE
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

    /**
     * @brief A
     * @param node
     * @param row
     * @param col
     * @return
     */
    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;

    /**
     * @brief B
     * @param node
     * @param row
     * @param col
     * @return
     */
    virtual auto B(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double = 0;

    /**
     * @brief C
     * @param node
     * @param row
     * @return
     */
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double = 0;

protected:
    /**
     * @brief initial
     * @param condition
     * @param row
     * @return
     */
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;

    /**
     * @brief boundary
     * @param node
     * @param condition
     * @param row
     * @return
     */
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double = 0;

protected:

    /**
     * @brief count
     * @return
     */
    virtual auto count() const -> unsigned int = 0;
};

/**
 * @brief The LinearBoundaryValueProblemODE class
 * r(x)y"(x)+p(x)y'(x)+q(x)y(x)=f(x);
 * y(0) = boundary(Left);
 * y(N) = boundart(Right);
 * @see BoundaryValueProblemODE
 * @see BoundaryValueProblem
 */
class MINIMUMSHARED_EXPORT LinearBoundaryValueProblemODE : protected BoundaryValueProblemODE
{
protected:
    virtual double r(unsigned int n) const = 0;
    virtual double p(unsigned int n) const = 0;
    virtual double q(unsigned int n) const = 0;
    virtual double f(unsigned int n) const = 0;

public:
    void calculateX(DoubleVector &x, double h, unsigned int N) const;

    void calculateN2L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN2R2LD(DoubleVector &x, double h, unsigned int N) const;

    void calculateN4L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN4R2LD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN4CNTR(DoubleVector &x, double h, unsigned int N) const;

    void calculateN6L2RD(DoubleVector &x, double h, unsigned int N) const;
    void calculateN6R2LD(DoubleVector &x, double h, unsigned int N) const;
};

#endif // SECOND_ORDER_LINEAR_ODE_H
