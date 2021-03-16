#ifndef TESTFUNCTION_H
#define TESTFUNCTION_H

#include "border_global.h"

#define FUNCTION_T1
#define FUNCTION_X1
#define FUNCTION_Y1

class BORDERSHARED_EXPORT TestFunction
{
public:
    enum Derivative
    {
        FunctionValue = 0,
        TimeFirstDerivative,
        TimeSecondDerivative,
        SpaceFirstDerivativeX,
        SpaceSecondDerivativeX,
        SpaceFirstDerivativeY,
        SpaceSecondDerivativeY,
        SpaceNorm
    };

    static double u(const TimeNodePDE &tn, const SpaceNodePDE &sn, Derivative derivative, const Dimension &dimX = 0, const Dimension &dimY = 0);
};

double BORDERSHARED_EXPORT test1();

#endif // TESTFUNCTION_H
