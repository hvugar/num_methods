#include "ibvp.h"

void InitialBoundaryValueProblemPDE::setGrid(const Grid &grid)
{
    mgrid = grid;
}

const Grid& InitialBoundaryValueProblemPDE::grid() const
{
    return mgrid;
}

