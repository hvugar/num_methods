#ifndef CAUCHYPROBLEMMEX_H
#define CAUCHYPROBLEMMEX_H

//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>
#include <printer.h>

#define SAMPLE_3

class CauchyProblemMEx : public NonLinearODE1stOrder
{
public:
    static void Main(int agrc, char *argv[]);
protected:
    virtual double f(double x, const DoubleVector &y, unsigned int k, unsigned int i) const;

    double y(double x, unsigned int k, unsigned int i) const;
};

#endif // CAUCHYPROBLEMMEX_H
