#include "ibvp.h"

InitialBoundaryValueProblemPDE::InitialBoundaryValueProblemPDE() {}

InitialBoundaryValueProblemPDE::InitialBoundaryValueProblemPDE(const InitialBoundaryValueProblemPDE& other) :
    _timeDimension(other._timeDimension), _spaceDimensionX(other._spaceDimensionX),
    _spaceDimensionY(other._spaceDimensionY), _spaceDimensionZ(other._spaceDimensionZ) {}

InitialBoundaryValueProblemPDE& InitialBoundaryValueProblemPDE::operator =(const InitialBoundaryValueProblemPDE& other)
{
    this->_timeDimension = other._timeDimension;
    this->_spaceDimensionX = other._spaceDimensionX;
    this->_spaceDimensionY = other._spaceDimensionY;
    this->_spaceDimensionZ = other._spaceDimensionZ;
    return *this;
}

InitialBoundaryValueProblemPDE::~InitialBoundaryValueProblemPDE() {}

auto InitialBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension) -> void
{
    _timeDimension = dimension;
}

auto InitialBoundaryValueProblemPDE::timeDimension() const -> const Dimension &
{
    return _timeDimension;
}

auto InitialBoundaryValueProblemPDE::setSpaceDimensionX(const Dimension &dimension) -> void
{
    _spaceDimensionX = dimension;
}

auto InitialBoundaryValueProblemPDE::spaceDimensionX() const -> const Dimension&
{
    return _spaceDimensionX;
}

auto InitialBoundaryValueProblemPDE::setSpaceDimensionY(const Dimension &dimension) -> void
{
    _spaceDimensionY = dimension;
}

auto InitialBoundaryValueProblemPDE::spaceDimensionY() const -> const Dimension&
{
    return _spaceDimensionY;
}

auto InitialBoundaryValueProblemPDE::setSpaceDimensionZ(const Dimension &dimension) -> void
{
    _spaceDimensionZ = dimension;
}

auto InitialBoundaryValueProblemPDE::spaceDimensionZ() const -> const Dimension&
{
    return _spaceDimensionZ;
}

auto InitialBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY) -> void
{
    this->_spaceDimensionX = dimensionX;
    this->_spaceDimensionY = dimensionY;
}

auto InitialBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension& dimensionZ) -> void
{
    this->_spaceDimensionX = dimensionX;
    this->_spaceDimensionY = dimensionY;
    this->_spaceDimensionZ = dimensionZ;
}
