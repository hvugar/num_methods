#ifndef CAUCHYPROBLEM_H
#define CAUCHYPROBLEM_H

#include "ibvp.h"

/**
 * @brief The CauchyProblem class
 * Ordinary differensial equation for of dy/dx=f(x,y).
 * Where x is independent variable. y(x) is searched function.
 */
class MINIMUMSHARED_EXPORT CauchyProblem : public InitialValueProblem
{
public:
    CauchyProblem(const Dimension &grid);
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

class MINIMUMSHARED_EXPORT CauchyProblemM : public InitialValueProblem
{
public:
    CauchyProblemM(const Dimension &grid);

    void calculate(double x0, const DoubleVector &y0, DoubleMatrix &ry, Method method = RK4, Direction direction = L2R);
    const Dimension &grid() const;

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const = 0;

private:
    void calculateRK2(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);

private:
    Dimension mgrid;
};

#endif // CAUCHYPROBLEM_H
