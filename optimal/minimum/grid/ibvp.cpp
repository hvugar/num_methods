#include "ibvp.h"
#include "linearequation.h"
#include <math.h>

InitialBoundaryValueProblemPDE::~InitialBoundaryValueProblemPDE() {}

auto InitialBoundaryValueProblemPDE::dimSize() const -> unsigned int
{
    return static_cast<unsigned>(mspaceDimension.size());
}

auto InitialBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension) -> void
{
    mtimeDimension = dimension;
}

auto InitialBoundaryValueProblemPDE::timeDimension() const -> const Dimension &
{
    return mtimeDimension;
}

auto InitialBoundaryValueProblemPDE::addSpaceDimension(const Dimension &dimension) -> void
{
    mspaceDimension.push_back(dimension);
}

auto InitialBoundaryValueProblemPDE::spaceDimension(Dimension::SpaceDimension dim) const -> const Dimension&
{
    return mspaceDimension.at(dim-1);
}
