#include "test_function.h"

#define FUNCTION_T1
#define FUNCTION_X3
#define FUNCTION_Y1

double TestFunction::u(const TimeNodePDE &tn, const SpaceNodePDE &sn, Derivative derivative, const Dimension &dimX, const Dimension &dimY)
{
    double res = 0.0;
    switch (derivative)
    {
    case FunctionValue:
    {
        res = 0.0;
#ifdef FUNCTION_T1
        res += tn.t;
#endif
#ifdef FUNCTION_T2
        res += tn.t*tn.t;
#endif
#ifdef FUNCTION_T3
        res += tn.t*tn.t*tn.t;
#endif
/******************************************/
#ifdef FUNCTION_X1
        res += sn.x;
#endif
#ifdef FUNCTION_X2
        res += sn.x*sn.x;
#endif
#ifdef FUNCTION_X3
        res += sn.x*sn.x*sn.x;
#endif
/******************************************/
#ifdef FUNCTION_Y1
        res += sn.y;
#endif
#ifdef FUNCTION_Y2
        res += sn.y*sn.y;
#endif
#ifdef FUNCTION_Y3
        res += sn.y*sn.y*sn.y;
#endif

    } break;
    case TimeFirstDerivative:
    {
#ifdef FUNCTION_T1
        res = 1.0;
#endif
#ifdef FUNCTION_T2
        res = 2.0*tn.t;
#endif
#ifdef FUNCTION_T3
        res = 3.0*tn.t*tn.t;
#endif
    } break;
    case TimeSecondDerivative:
    {
#ifdef FUNCTION_T1
        res = 0.0;
#endif
#ifdef FUNCTION_T2
        res = 2.0;
#endif
#ifdef FUNCTION_T3
        res = 6.0*tn.t;
#endif
    } break;
    case SpaceFirstDerivativeX:
    {
#ifdef FUNCTION_X1
        res = 1.0;
#endif
#ifdef FUNCTION_X2
        res = 2.0*sn.x;
#endif
#ifdef FUNCTION_X3
        res = 3.0*sn.x*sn.x;
#endif
    } break;
    case SpaceSecondDerivativeX:
    {
#ifdef FUNCTION_X1
        res = 0.0;
#endif
#ifdef FUNCTION_X2
        res = 2.0;
#endif
#ifdef FUNCTION_X3
        res = 6.0*sn.x;
#endif
    } break;
    case SpaceFirstDerivativeY:
    {
#ifdef FUNCTION_Y1
        res = 1.0;
#endif
#ifdef FUNCTION_Y2
        res = 2.0*sn.y;
#endif
#ifdef FUNCTION_Y3
        res = 3.0*sn.y*sn.y;
#endif
    } break;
    case SpaceSecondDerivativeY:
    {
#ifdef FUNCTION_Y1
        res = 0.0;
#endif
#ifdef FUNCTION_Y2
        res = 2.0;
#endif
#ifdef FUNCTION_y3
        res = 6.0*sn.y;
#endif
    } break;
    case SpaceNorm:
    {
#ifdef FUNCTION_X1
        if (sn.i == dimX.min()) res = -1.0;
        if (sn.i == dimX.max()) res = +1.0;
#endif
#ifdef FUNCTION_X2
        if (sn.i == dimX.min()) res = -2.0*sn.x;
        if (sn.i == dimX.max()) res = +2.0*sn.x;
#endif
#ifdef FUNCTION_X3
        if (sn.i == dimX.min()) res = -3.0*sn.x*sn.x;
        if (sn.i == dimX.max()) res = +3.0*sn.x*sn.x;
#endif
#ifdef FUNCTION_Y1
        if (sn.j == dimY.min()) res = -1.0;
        if (sn.j == dimY.max()) res = +1.0;
#endif
#ifdef FUNCTION_Y2
        if (sn.j == dimY.min()) res = -2.0*sn.y;
        if (sn.j == dimY.max()) res = +2.0*sn.y;
#endif
#ifdef FUNCTION_Y3
        if (sn.j == dimY.min()) res = -3.0*sn.y*sn.y;
        if (sn.j == dimY.max()) res = +3.0*sn.y*sn.y;
#endif
    } break;
    }

    return res;
}
