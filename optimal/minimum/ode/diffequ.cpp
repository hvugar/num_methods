#include "diffequ.h"

DifferentialEquation::DifferentialEquation() {}

DifferentialEquation::~DifferentialEquation() {}

OrdinaryDifferentialEquation::OrdinaryDifferentialEquation() {}

OrdinaryDifferentialEquation::~OrdinaryDifferentialEquation() {}

unsigned int DifferentialEquation::equationsNumber() const
{
    return 1;
}

const Dimension& OrdinaryDifferentialEquation::dimension() const
{
    return _dimension;
}

void OrdinaryDifferentialEquation::setDimension(const Dimension &dimension)
{
    _dimension = dimension;
}

const Dimension& SystemDifferentialEquation::dimension() const
{
    return _dimension;
}

void SystemDifferentialEquation::setDimension(const Dimension &dimension)
{
    _dimension = dimension;
}
