#include "diffequ.h"

unsigned int DifferentialEquation::equationsNumber() const
{
    return 1;
}

const ODEGrid& DifferentialEquation::grid() const
{
    return mgrid;
}

void DifferentialEquation::setGrid(const ODEGrid& grid)
{
    mgrid = grid;
}

const ODEGrid& SystemDifferentialEquation::grid() const
{
    return mgrid;
}

void SystemDifferentialEquation::setGrid(const ODEGrid &grid)
{
    mgrid = grid;
}
