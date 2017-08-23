#ifndef NONLINEARODE1STORDER_H
#define NONLINEARODE1STORDER_H

#include "diffequ.h"

/**
 * @brief The NonLinear ODE1 1st order in canonical (normal) form y'(x) = f(x, y(x));
 */
class MINIMUMSHARED_EXPORT NonLinearODE1stOrder : virtual public NonLinearODE, virtual public ODE1stOrder
{
public:
    void cauchyProblem(double x0, double y0, DoubleVector &ry, Method method = RK4, Direction direction = L2R);
public:
    void cauchyProblem(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Method method = RK4, Direction direction = L2R);
public:
    void cauchyProblem(double x0, const DoubleVector &y0, DoubleVector &ry, Method method = RK4, Direction direction = L2R);

protected:
    virtual double f(double x, double y, unsigned int k) const;

    virtual double f(double x, DoubleVector y, unsigned int k, unsigned int i) const;

private:
    void calculateRK2(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateRK4(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateEuler(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateEulerMod(double x0, double y0, DoubleVector &ry, Direction direction = L2R);

private:
    void calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
private:
    void calculateRK4(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
};
#endif // NONLINEARODE1STORDER_H
