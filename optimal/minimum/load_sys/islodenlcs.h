#ifndef ISYSTEMLINEARODENONLOCALCONTIONS_H
#define ISYSTEMLINEARODENONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>
#include <vector>
//#include <ode/cauchyp.h>
#include <ode/nlode1o.h>

class MINIMUMSHARED_EXPORT ISystemLinearODENonLocalContions
{
public:
    enum ConditionType
    {
        SeparatedLeft = 0,
        SeparatedRight = 1,
        NonSeparated = 2
    };

    struct Condition
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix mtrx;
        unsigned int index;
    };

    struct LoadPoint
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix mtrx;
    };
};

#endif // ISYSTEMLINEARODENONLOCALCONTIONS_H
