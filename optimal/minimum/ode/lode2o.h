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
    void solveBoundaryValueProblem(const DoubleVector &left, const DoubleVector &right, std::vector<DoubleVector> &ry);
    /**
     * @brief solveCauchyProblem
     * @param initial Initial conditions
     * @param initder Initial derivative conditions
     * @param ry      Result functions
     */
    void solveCauchyProblem(const DoubleVector &initial, const DoubleVector &initder, std::vector<DoubleVector> &ry);

protected:
    /**
     * @brief P one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double R(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief P one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double P(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief Q one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double Q(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief R one dimensional vector-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double F(const PointNodeODE &node, unsigned int row = 0) const = 0;

private:

};

#endif // SECOND_ORDER_LINEAR_ODE_H
