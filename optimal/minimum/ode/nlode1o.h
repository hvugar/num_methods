#ifndef NON_LINEAR_ODE1ST_ORDER_H
#define NON_LINEAR_ODE1ST_ORDER_H

#include "diffequ.h"

/**
 * @brief The NonLinear ODE1 1st order in canonical (normal) form y'(x) = f(x, y(x));
 */
class MINIMUMSHARED_EXPORT NonLinearODE1stOrder : virtual public NonLinearODE
{
public:
    void cauchyProblem(double x0, double y0, DoubleVector &ry, Method method = RK4, Direction direction = L2R);
    void cauchyProblem(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Method method = RK4, Direction direction = L2R);
    void cauchyProblem(double x0, const DoubleVector &y0, DoubleVector &ry, Method method = RK4, Direction direction = L2R);

protected:
    virtual double f(double x, double y, unsigned int k) const;
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;

private:
    void calculateRK2(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateRK4(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateEuler(double x0, double y0, DoubleVector &ry, Direction direction = L2R);
    void calculateEulerMod(double x0, double y0, DoubleVector &ry, Direction direction = L2R);

    void calculateRK2(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R);

    void calculateRK2(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
    void calculateRK4(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
    void calculateEuler(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
    void calculateEulerMod(double x0, const DoubleVector &y0, DoubleVector &y, Direction direction = L2R);
};
#endif // NON_LINEAR_ODE1ST_ORDER_H
