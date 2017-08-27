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
        DoubleMatrix mtrx;
        unsigned int nmbr;
        unsigned int index;
    };

    void calculate(const std::vector<Condition> &nscs, const DoubleVector &bt, std::vector<DoubleVector> &x);

    /* high order accuracy */

//private:
    void highOder2Accuracy(const std::vector<Condition>& cs, const DoubleVector& rs, std::vector<DoubleVector>& x);
    void highOder4Accuracy(const std::vector<Condition>& cs, const DoubleVector& rs, std::vector<DoubleVector>& x);
    void highOder6Accuracy(const std::vector<Condition>& cs, const DoubleVector& rs, std::vector<DoubleVector>& x);
    void discretisation(const std::vector<Condition>& cs, double* b) const;

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
