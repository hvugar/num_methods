#ifndef INITIALBOUNDARYVALUEPROBLEM_H
#define INITIALBOUNDARYVALUEPROBLEM_H

#include "bvp.h"
#include "grid.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : protected BoundaryValueProblemPDE, protected InitialValueProblem
{
public:
//    void setGrid(const GridPDE &grid);
//    const GridPDE& grid() const;

protected:
//    GridPDE mgrid;
};

#endif // INITIALBOUNDARYVALUEPROBLEM_H
