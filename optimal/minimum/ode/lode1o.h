#ifndef LINEAR_ODE1ST_ORDER_H
#define LINEAR_ODE1ST_ORDER_H

#include "diffequ.h"
#include "../matrix2d.h"
#include <printer.h>

/**
 * @brief The Linear ODE 1st order in canonical (normal) form y'(x) = A(x)y(x) + B(x);
 */
class MINIMUMSHARED_EXPORT LinearODE1stOrder : virtual public LinearODE, virtual public ODE1stOrder
{
public:
    struct Condition
    {
        double time;
        unsigned int nmbr;
        DoubleMatrix mtrx;
        unsigned int index;
    };

    void calculate(const std::vector<Condition> &nscs, const DoubleVector &bt, std::vector<DoubleVector> &x);

    void highOder2Accuracy(const std::vector<Condition> &cnd, const DoubleVector& rs);
    void highOder4Accuracy();
    void highOder6Accuracy();

protected:
    /**
     * @brief A one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double A(double x, unsigned int i, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief B one dimensional vector-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double B(double x, unsigned int i, unsigned int row = 0) const = 0;
};


#endif // LINEAR_ODE1ST_ORDER_H
