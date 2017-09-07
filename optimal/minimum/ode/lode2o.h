#ifndef LINEARODE2NDORDER_H
#define LINEARODE2NDORDER_H

#include "diffequ.h"

/**
 * @brief The Linear ODE 2nd order in canonical (normal) form y"(x) = P(x)y'(x) + Q(x)y(x)+R(x);
 */
class MINIMUMSHARED_EXPORT LinearODE2ndOrder : virtual public LinearODE
{
protected:
    /**
     * @brief P one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double P(double x, unsigned int i, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief Q one dimensional matrix-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double Q(double x, unsigned int i, unsigned int row = 0, unsigned int col = 0) const = 0;
    /**
     * @brief R one dimensional vector-function
     * @param x independent variable
     * @param i index of independent variable of given grid
     * @return
     */
    virtual double R(double x, unsigned int i, unsigned int row = 0) const = 0;
};

#endif // LINEARODE2NDORDER_H
