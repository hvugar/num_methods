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

#endif // RUNGEKUTTA_H
