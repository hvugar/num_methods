#include "ibvp.h"
#include "linearequation.h"

InitialBoundaryValueProblemPDE::~InitialBoundaryValueProblemPDE() {}

unsigned int InitialBoundaryValueProblemPDE::dimSize()
{
    return mspaceDimension.size();
}

void InitialBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension)
{
    mtimeDimension = dimension;
}

const Dimension &InitialBoundaryValueProblemPDE::timeDimension() const
{
    return mtimeDimension;
}

void InitialBoundaryValueProblemPDE::addSpaceDimension(const Dimension &dimension)
{
    mspaceDimension.push_back(dimension);
}

const Dimension &InitialBoundaryValueProblemPDE::spaceDimension(Dimension::SpaceDimension dim) const
{
    return mspaceDimension.at(dim);
}
