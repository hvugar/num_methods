#ifndef CAUCHYPROBLEM_H
#define CAUCHYPROBLEM_H

#include "ibvp.h"

class MINIMUMSHARED_EXPORT CauchyProblem1 : public InitialValueProblem
{
public:
    CauchyProblem1(const Dimension &grid);
    void calculate(double x0, double y0, DoubleVector &y, Method method = RK4, Direction direction = L2R);

protected:
    virtual double f(double xn, double yn) const = 0;

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
    void calculate(double x0, const DoubleVector &y0, DoubleMatrix &y, Method method = RK4, Direction direction = L2R);

protected:
    virtual double f(double x, const DoubleVector &y, unsigned int i) const = 0;

private:
    void calculateRK2(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, DoubleMatrix &y, Direction direction = L2R);

private:
    Dimension mgrid;
};

#endif // CAUCHYPROBLEM_H
