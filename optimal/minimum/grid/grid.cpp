#include "grid.h"

Dimension::Dimension(double step, double min, double max, unsigned int minN, unsigned int maxN)
    : mstep(step), mmin(min), mmax(max), mminN(minN), mmaxN(maxN)
{}

double Dimension::step() const { return mstep; }

double Dimension::min() const { return mmin; }

double Dimension::max() const { return mmax; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

SpaceDimension::SpaceDimension(double step, double min, double max, unsigned int minN, unsigned int maxN)
    : Dimension(step, min, max, minN, maxN)
{}

TimeDimension::TimeDimension(double step, double min, double max, unsigned int minN, unsigned int maxN)
    : Dimension(step, min, max, minN, maxN)
{}

const TimeDimension& GridPDE::timeDimension() const
{
    return mTimeDimension;
}

void GridPDE::setTimeDimension(const TimeDimension &timeDimension)
{
    mTimeDimension = timeDimension;
}

const SpaceDimension& GridPDE::spaceDimensions(SpaceDimension::DimensionSize dim) const
{
    return mSpaceDimensions.at(dim);
}

void GridPDE::addSpaceDimension(const SpaceDimension &spaceDimension)
{
    mSpaceDimensions.push_back(spaceDimension);
}

unsigned int GridPDE::spaceDimSize() const
{
    return mSpaceDimensions.size();
}

const Dimension& GridODE::dimension() const
{
    return mdimension;
}

void GridODE::setDimension(const Dimension &dimension)
{
    mdimension = dimension;
}
