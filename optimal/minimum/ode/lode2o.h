#ifndef SECOND_ORDER_LINEAR_ODE_H
#define SECOND_ORDER_LINEAR_ODE_H

#include "diffequ.h"

/**
 * @brief Линейное дифференциальное уравнение второго порядка с переменными коэффициентами
 * The Linear ODE 2nd order in canonical (normal) form R(x)y"(x) + P(x)y'(x) + Q(x)y(x) = F(x);
 */
class MINIMUMSHARED_EXPORT SecondOrderLinearODE : virtual public LinearODE
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
     * @brief solveCauchyProblem
     * @param initial Initial conditions
     * @param initder Initial derivative conditions
     * @param ry      Result functions
     */
    void solveCauchyProblem(const DoubleVector &initial, const DoubleVector &initder, std::vector<DoubleVector> &ry) const;

protected:
    virtual double A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const = 0;
    virtual double B(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const = 0;
    virtual double C(const PointNodeODE &node, unsigned int row = 1) const = 0;

private:

};

#endif // SECOND_ORDER_LINEAR_ODE_H
