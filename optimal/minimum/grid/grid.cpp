#include "grid.h"
#include <cstdio>

UniformODEGrid::UniformODEGrid(double step, int min, int max)
    : mstep(step), mminN(min), mmaxN(max)
{}

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

/////////////////////////////////////////////////////////////////////////////////////////////////////

UniformPDEGrid::UniformPDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaceDimensions UNUSED_PARAM)
{
    mtimeDimension = timeDimension;
    mspaceDimensions = spaceDimensions;
}

const Dimension &UniformPDEGrid::timeDimension() const
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
