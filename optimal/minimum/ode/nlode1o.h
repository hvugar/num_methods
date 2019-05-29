#ifndef FIRTS_ORDER_NON_LINEAR_ODE_H
#define FIRTS_ORDER_NON_LINEAR_ODE_H

#include "diffequ.h"
#include "../grid/ibvp.h"

/**
 * @brief The FirstOrderNonLinearODE in canonical (normal) form y'(x) = f(x, y(x));
 */
class MINIMUMSHARED_EXPORT FirstOrderNonLinearODE : virtual public NonLinearODE
{
public:
    void solveInitialValueProblem(DoubleVector &rv, ODESolverMethod method = ODESolverMethod::EULER) const;
    void solveInitialValueProblem(std::vector<DoubleVector> &rv, ODESolverMethod method = ODESolverMethod::EULER) const;

    void cauchyProblem(double x0, double y0, DoubleVector &ry, ODESolverMethod method = ODESolverMethod::RUNGE_KUTTA_4, Direction direction = Direction::L2R);
    void cauchyProblem(double x0, const DoubleVector &y0, DoubleVector &ry, ODESolverMethod method = ODESolverMethod::RUNGE_KUTTA_4, Direction direction = Direction::L2R);
    void cauchyProblem(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, ODESolverMethod method = ODESolverMethod::RUNGE_KUTTA_4, Direction direction = Direction::L2R);

protected:
    virtual double f(double x, double y, unsigned int k) const;
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;

    /**
     * @brief f
     * @param node
     * @param x
     * @param m
     * @return
     */
    virtual double f(const PointNodeODE &node, const DoubleVector &x, unsigned int m = 1) const = 0;

protected:
    //virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double = 0;
    //virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double = 0;

private:
    void calculateRK2(double x0, double y0, DoubleVector &ry, Direction direction = Direction::L2R);
    void calculateRK4(double x0, double y0, DoubleVector &ry, Direction direction = Direction::L2R);
    void calculateEuler(double x0, double y0, DoubleVector &ry, Direction direction = Direction::L2R);
    void calculateEulerMod(double x0, double y0, DoubleVector &ry, Direction direction = Direction::L2R);

    void calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = Direction::L2R);
    void calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = Direction::L2R);
    void calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = Direction::L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = Direction::L2R);

    void calculateRK2(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = Direction::L2R);
    void calculateRK4(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = Direction::L2R);
    void calculateEuler(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = Direction::L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = Direction::L2R);

    auto solveInitialValueProblemEuler(DoubleVector &rv) const -> void;
    auto solveInitialValueProblemEulerMod(DoubleVector &rv) const -> void;
    auto solveInitialValueProblemRK2(DoubleVector &rv) const -> void;
    auto solveInitialValueProblemRK4(DoubleVector &rv) const -> void;
    auto solveInitialValueProblemRK6(DoubleVector &rv) const -> void;

    auto solveInitialValueProblemEuler(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemEulerMod(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK2(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK4(std::vector<DoubleVector> &rv) const -> void;
    auto solveInitialValueProblemRK6(std::vector<DoubleVector> &rv) const -> void;

};

#endif // FIRTS_ORDER_NON_LINEAR_ODE_H
