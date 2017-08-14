#include "uniformgrid.h"

UniformGrid::UniformGrid(double step, int min, int max) :
    mstep(step), mminN(min), mmaxN(max) {}

double UniformGrid::step() const
{
    return mstep;
}

int UniformGrid::minN() const
{
    return mminN;
}

int UniformGrid::maxN() const
{
    return mmaxN;
}

int UniformGrid::sizeN() const
{
    return mmaxN - mminN;
}

bool UniformGrid::isGridSet() const
{
    return mminN != 0 || mmaxN != 0;
}

