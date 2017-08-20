#include "diffequ.h"

const ODEGrid& SystemDifferentialEquation::grid() const
{
    return mgrid;
}

void SystemDifferentialEquation::setGrid(const ODEGrid &grid)
{
    mgrid = grid;
}
