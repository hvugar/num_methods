#include "islodenlcs.h"

void ISystemLinearODENonLocalContions::setGrid(const ODEGrid &grid)
{
    this->mgrid = grid;
}

const ODEGrid& ISystemLinearODENonLocalContions::grid() const
{
    return mgrid;
}

