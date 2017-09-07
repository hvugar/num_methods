#include "diffequ.h"

unsigned int DifferentialEquation::equationsNumber() const
{
    return 1;
}

const UniformODEGrid& DifferentialEquation::grid() const
{
    return mgrid;
}

void DifferentialEquation::setGrid(const UniformODEGrid& grid)
{
    mgrid = grid;
}

const UniformODEGrid& SystemDifferentialEquation::grid() const
{
    return mgrid;
}

void SystemDifferentialEquation::setGrid(const UniformODEGrid &grid)
{
    mgrid = grid;
}
