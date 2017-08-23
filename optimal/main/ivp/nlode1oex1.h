#ifndef NON_LINEAR_ODE1ST_ORDER_EX1_H
#define NON_LINEAR_ODE1ST_ORDER_EX1_H

//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>
#include <printer.h>

#define SAMPLE_3

class NonLinearODE1stOrderEx1 : public NonLinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);
protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;

    double y(double x, unsigned int k, unsigned int i) const;
};

#endif // NON_LINEAR_ODE1ST_ORDER_EX1_H
