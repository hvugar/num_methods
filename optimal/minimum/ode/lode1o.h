#ifndef FIRTS_ORDER_LINEAR_ODE_H
#define FIRTS_ORDER_LINEAR_ODE_H

#include "diffequ.h"

struct MINIMUMSHARED_EXPORT NonLocalCondition
{
    NonLocalCondition();
    NonLocalCondition(unsigned int i, const PointNodeODE &node, const DoubleMatrix &m);
    virtual ~NonLocalCondition();

    unsigned int i;
    PointNodeODE n;
    DoubleMatrix m;
};

struct Condition
{
    double time;
    DoubleMatrix mtrx;
    unsigned int nmbr;
    unsigned int index;
};

/**
 * @brief Линейное дифференциальное уравнение первого порядка с переменными коэффициентами
 * The Linear ODE 1st order in canonical (normal) form y'(x) = A(x)y(x) + B(x);
 */
class MINIMUMSHARED_EXPORT FirstOrderLinearODE : virtual public LinearODE
{
public:
    enum class AccuracyStep
    {
        Step_2 = 2,
        Step_4 = 4,
        Step_6 = 6
    };

    /**
     * @brief solve transfer of conditions
     * @param C
     * @param d
     * @param x
     * @param k
     * @param M
     * @param direction
     */
    void transferOfCondition(const std::vector<NonLocalCondition> &C, const DoubleVector &d, std::vector<DoubleVector> &x, unsigned int k) const;

protected:
    /**
     * @brief A  A nxn dimensional matrix-function
     * @param node
     * @param row <= n
     * @param col <= n
     * @return
     */

    virtual double A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const = 0;
    /**
     * @brief B n dimensional vector-function
     * @param node
     * @param row
     * @return
     */
    virtual double B(const PointNodeODE &node, unsigned int row = 1) const = 0;

private:
    void discritize(const std::vector<NonLocalCondition> &co, std::vector<NonLocalCondition> &cn, unsigned int k=4) const;
};


#endif // FIRTS_ORDER_LINEAR_ODE_H
