#include "ibvp.h"

void InitialBoundaryValueProblemPDE::setGrid(const GridPDE &grid)
{
    mgrid = grid;
}

const GridPDE& InitialBoundaryValueProblemPDE::grid() const
{
    return mgrid;
}

