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
        SpaceFirstDerivativeX1,
        SpaceSecondDerivativeX1,
        SpaceFirstDerivativeX2,
        SpaceSecondDerivativeX2,
    };

    static double u(const TimeNodePDE &tn, const SpaceNodePDE &sn, Derivative derivative);
};

#endif // TESTFUNCTION_H
