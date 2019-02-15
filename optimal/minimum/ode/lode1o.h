#ifndef LINEAR_ODE1ST_ORDER_H
#define LINEAR_ODE1ST_ORDER_H

#include "diffequ.h"

struct NonLocalCondition
{
    double x;
    DoubleMatrix m;
    unsigned int i;
};

/**
 * @brief Линейное дифференциальное уравнение первого порядка с переменными коэффициентами
 * The Linear ODE 1st order in canonical (normal) form y'(x) = A(x)y(x) + B(x);
 */
class MINIMUMSHARED_EXPORT LinearODE1stOrder : virtual public LinearODE
{
public:
    struct Condition
    {
        double time;
        DoubleMatrix mtrx;
        unsigned int nmbr;
        unsigned int index;
    };

    void calculate(const std::vector<Condition> &cs, const DoubleVector &bt, std::vector<DoubleVector> &x);

    void solveHighOderAccuracy(const std::vector<Condition>& cs, const DoubleVector& rs, std::vector<DoubleVector>& x, unsigned int k, Direction direction = L2R);

    void calculate(double x0, const DoubleVector &y0, std::vector<DoubleVector> &ry, Direction direction = L2R) const;
    void calculate(double x0, double y0, std::vector<double> &ry, Direction direction = L2R) const;

    void solve(const std::vector<NonLocalCondition> &C, const DoubleVector &d,
               std::vector<DoubleVector> &x,
               unsigned int k, unsigned int M, Direction direction = Direction::L2R) const;

    /* high order accuracy */

private:
    void highOder2Accuracy(const std::vector<Condition> &cs, const DoubleVector &rs, std::vector<DoubleVector> &x, Direction direction = L2R);
    void highOder4Accuracy(const std::vector<Condition> &cs, const DoubleVector &rs, std::vector<DoubleVector> &x, Direction direction = L2R);
    void highOder6Accuracy(const std::vector<Condition> &cs, const DoubleVector &rs, std::vector<DoubleVector> &x, Direction direction = L2R);

protected:
    /**
     * @brief A one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double A(const PointNodeODE &node, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief B one dimensional vector-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double B(const PointNodeODE &node, unsigned int row = 0) const = 0;
};


#endif // LINEAR_ODE1ST_ORDER_H
