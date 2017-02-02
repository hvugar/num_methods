#include "grid.h"

Dimension::Dimension(double step, unsigned int maxN, unsigned int minN)
    : mstep(step), mmaxN(maxN), mminN(minN)
{
    mmin = mminN*mstep;
    mmax = mmaxN*mstep;
}

double Dimension::step() const { return mstep; }

double Dimension::min() const { return mmin; }

double Dimension::max() const { return mmax; }

unsigned int Dimension::minN() const { return mminN; }

unsigned int Dimension::maxN() const { return mmaxN; }

unsigned int Dimension::sizeN() const { return mmaxN-mminN; }

SpaceDimension::SpaceDimension(double step, unsigned int maxN, unsigned int minN) : Dimension(step, maxN, minN)
{}

TimeDimension::TimeDimension(double step, unsigned int maxN, unsigned int minN) : Dimension(step, maxN, minN)
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
