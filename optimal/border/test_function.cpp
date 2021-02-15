#include "test_function.h"

double TestFunction::u(const TimeNodePDE &tn, const SpaceNodePDE &sn, Derivative derivative)
{
    double res = 0.0;
    switch (derivative)
    {
    case FunctionValue:
    {
        res = sn.x + sn.y + tn.t;
    }
        break;
    case TimeFirstDerivative: res = 1.0; break;
    case TimeSecondDerivative: res = 0.0; break;
    case SpaceFirstDerivativeX1: res = 1.0; break;
    case SpaceSecondDerivativeX1: res = 0.0; break;
    case SpaceFirstDerivativeX2: res = 1.0; break;
    case SpaceSecondDerivativeX2: res = 0.0; break;
    }

    return res;
}
