#include "grid.h"
#include <cstdio>

ODEGrid::ODEGrid() {}

ODEGrid::ODEGrid(const Dimension &dimension) : mdimension(dimension) {}

const Dimension& ODEGrid::dimension() const
{
    return mdimension;
}

void ODEGrid::setDimension(const Dimension &dimension)
{
    mdimension = dimension;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

PDEGrid::PDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaces UNUSED_PARAM)
{
    mtimeDimension = timeDimension;
}

const Dimension &PDEGrid::timeDimension() const
{
    return mtimeDimension;
}

Dimension::Dimension(double step, unsigned int maxN, unsigned int minN)
    : mstep(step), mmaxN(maxN), mminN(minN)
{}

double Dimension::step() const { return mstep; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

unsigned int Dimension::sizeN() const { return mmaxN-mminN; }
