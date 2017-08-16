#ifndef CAUCHYPROBLEM_H
#define CAUCHYPROBLEM_H

#include "../grid/ibvp.h"
#include "diffequ.h"

/**
 * @brief The CauchyProblem class
 * Ordinary differensial equation for of dy/dx=f(x,y).
 * Where x is independent variable. y(x) is searched function.
 */
class MINIMUMSHARED_EXPORT CauchyProblem1stOrder : public InitialValueProblem
{
public:
    CauchyProblem1stOrder(const Dimension &grid);
    /**
     * @brief calculate
     * @param x0 initial value for independent variable x.
     * @param y0 initial value of searched y(x) function at point x0. y(x0) = y0.
     * @param y  Vector of discrete values of searched function in given grid.
     * @param method Using method for calculation of function values y(x).
     * @param direction Direction for finding function values which indicates that initial value are given start of end of interval.
     */
    void calculate(double x0, double y0, DoubleVector &y, Method method = RK4, Direction direction = L2R);

protected:
    /**
     * @brief f Right side of Cauchy problem dy/dx=f(x,y). Must be override.
     * @param x independend variable.
     * @param y searched function.
     * @param k index of grid.
     */
    virtual double f(double x, double y, unsigned int k) const = 0;

private:
    void calculateRK2(double x0, double y0, DoubleVector &y, Direction direction = L2R);
    void calculateRK4(double x0, double y0, DoubleVector &y, Direction direction = L2R);
    void calculateEuler(double x0, double y0, DoubleVector &y, Direction direction = L2R);
    void calculateEulerMod(double x0, double y0, DoubleVector &y, Direction direction = L2R);

private:
    Dimension mgrid;
};

class MINIMUMSHARED_EXPORT CauchyProblemM1stOrder : public SystemNonLinearODE1stOrder, public InitialValueProblem
{
public:
    CauchyProblemM1stOrder(const ODEGrid &grid);

public:
    void calculateCP(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Method method = RK4, Direction direction = L2R);
private:
    void calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);

public:
    void calculateCP(double x0, const DoubleVector &y0, DoubleVector &ry, Method method = RK4, Direction direction = L2R);
private:
    void calculateRK4(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
};

#endif // CAUCHYPROBLEM_H
