#ifndef INITIALBOUNDARYVALUEPROBLEM_H
#define INITIALBOUNDARYVALUEPROBLEM_H

#include "bvp.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : protected BoundaryValueProblemPDE
{
public:
    void setGrid(const Grid &grid);
    const Grid& grid() const;

private:
    Grid mgrid;
    BoundaryCondition condition;
};

#endif // INITIALBOUNDARYVALUEPROBLEM_H
