#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "global.h"
#include "function.h"

typedef double (*R2FunctionX)(double x, double y);

class MINIMUMSHARED_EXPORT RungeKutta
{
public:
    static void calculate(R2Function *f, double x0, double x1, double y0, double &y1, double dx);
    static unsigned int calculate(R2Function *f, double x0, double x1, double y0, DoubleVector& y, double dx);
    static void calculate(R2Function *f, double x0, double y0, DoubleVector& y, double dx);

    static void calculate(R2FunctionX f, double x0, double y0, DoubleVector& y, double dx);
};

class MINIMUMSHARED_EXPORT CauchyProblem : public R2Function
{
    virtual double fx(double x, const DoubleVector &y) const = 0;

    static void rungeKutta(const CauchyProblem *cp, double x0, double y0, double h, unsigned int N, DoubleVector &y);

private:
    double mx0;
    double my0;
};

#endif // RUNGEKUTTA_H
