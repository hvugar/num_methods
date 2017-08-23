#include "uniformgrid.h"

UniformODEGrid::UniformODEGrid(double step, int min, int max) : mstep(step), mminN(min), mmaxN(max) {}

double UniformODEGrid::step() const
{
    return mstep;
}

int UniformODEGrid::minN() const
{
    return mminN;
}

int UniformODEGrid::maxN() const
{
    return mmaxN;
}

int UniformODEGrid::sizeN() const
{
    return mmaxN - mminN;
}

bool UniformODEGrid::isGridSet() const
{
    return mminN != 0 || mmaxN != 0;
}

