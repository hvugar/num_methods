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

//auto InitialBoundaryValueProblemPDE::dimSize() const -> unsigned int
//{
//    return static_cast<unsigned>(mspaceDimension.size());
//}

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

//auto InitialBoundaryValueProblemPDE::addSpaceDimension(const Dimension &dimension) -> void
//{
//    mspaceDimension.push_back(dimension);
//}

//auto InitialBoundaryValueProblemPDE::spaceDimension(Dimension::SpaceDimension dim) const -> const Dimension&
//{
//    return mspaceDimension.at(dim-1);
//}
