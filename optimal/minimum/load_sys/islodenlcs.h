#ifndef ISYSTEMLINEARODENONLOCALCONTIONS_H
#define ISYSTEMLINEARODENONLOCALCONTIONS_H

#include <vector2d.h>
#include <matrix2d.h>
#include <grid/grid.h>
#include <vector>
#include <ode/cauchyp.h>

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
        std::vector<std::vector<DoubleVector>> vectors;
    };

    struct LoadPoint
    {
        ConditionType type;
        double time;
        unsigned int nmbr;
        DoubleMatrix mtrx;
    };

    void setGrid(const ODEGrid &grid);
    const ODEGrid& grid() const;

private:
    ODEGrid mgrid;
};

#endif // ISYSTEMLINEARODENONLOCALCONTIONS_H
