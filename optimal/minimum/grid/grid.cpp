#include "grid.h"
#include <cstdio>

UniformODEGrid::UniformODEGrid(double step, int min, int max) : mstep(step), mminN(min), mmaxN(max)
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

UniformPDEGrid::UniformPDEGrid()
{}

UniformPDEGrid::UniformPDEGrid(const Dimension &timeDimension, std::vector<Dimension> &spaceDimensions UNUSED_PARAM)
{
    mtimeDimension = timeDimension;
    mspaceDimensions = spaceDimensions;
}

void UniformPDEGrid::setTimeDimension(const Dimension &timeDimension)
{
    mtimeDimension = timeDimension;
}

const Dimension& UniformPDEGrid::timeDimension() const
{
    return mtimeDimension;
}

void UniformPDEGrid::addSpaceDimension(const Dimension &spaceDimension)
{
    mspaceDimensions.push_back(spaceDimension);
}

const Dimension& UniformPDEGrid::spaceDimension(Dimension::SpaceDimension dimension) const
{
    return mspaceDimensions[dimension];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

Dimension::Dimension(double step, int minN, int maxN) : mstep(step), mminN(minN), mmaxN(maxN)
{}

double Dimension::step() const { return mstep; }

int Dimension::minN() const { return mminN; }

int Dimension::maxN() const { return mmaxN; }

int Dimension::sizeN() const { return mmaxN-mminN; }
