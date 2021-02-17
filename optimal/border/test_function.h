#ifndef TESTFUNCTION_H
#define TESTFUNCTION_H

#include "border_global.h"

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

#endif // TESTFUNCTION_H
