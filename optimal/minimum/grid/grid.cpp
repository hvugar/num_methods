#include "grid.h"
#include <cstdio>

SpacePoint::SpacePoint(double x, double y, double z) : x(x), y(y), z(z) {}

SpaceNodePDE::SpaceNodePDE(int i, double x, int j, double y, int k, double z) : i(i), j(j), k(k)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

SpaceNodePDE::SpaceNodePDE(const SpaceNodePDE &sn) : i(sn.i), j(sn.j), k(sn.k)
{
    this->x = sn.x;
    this->y = sn.y;
    this->z = sn.z;
}

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

Dimension::Dimension(double step, int minN, int maxN) : m_step(step), m_min(minN), m_max(maxN)
{}

double Dimension::step() const { return m_step; }

int Dimension::min() const { return m_min; }

int Dimension::max() const { return m_max; }

int Dimension::size() const { return m_max-m_min; }
