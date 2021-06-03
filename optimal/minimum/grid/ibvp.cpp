#include "ibvp.h"

GridDimensionPDE::GridDimensionPDE() {}

GridDimensionPDE::GridDimensionPDE(const GridDimensionPDE& /*other*/) {}

GridDimensionPDE& GridDimensionPDE::operator=(const GridDimensionPDE& /*other*/) { return *this; }

GridDimensionPDE::~GridDimensionPDE() {}


InitialBoundaryValueProblemPDE::InitialBoundaryValueProblemPDE() {}

InitialBoundaryValueProblemPDE::InitialBoundaryValueProblemPDE(const InitialBoundaryValueProblemPDE& other) :
    InitialValueProblemPDE(other)/*, BoundaryValueProblemPDE(other),
      _timeDimension(other._timeDimension), _spaceDimensionX(other._spaceDimensionX),
      _spaceDimensionY(other._spaceDimensionY), _spaceDimensionZ(other._spaceDimensionZ)*/ {}

InitialBoundaryValueProblemPDE& InitialBoundaryValueProblemPDE::operator =(const InitialBoundaryValueProblemPDE& other)
{
    C_UNUSED(other);
    //this->_timeDimension = other._timeDimension;
    //this->_spaceDimensionX = other._spaceDimensionX;
    //this->_spaceDimensionY = other._spaceDimensionY;
    //this->_spaceDimensionZ = other._spaceDimensionZ;
    return *this;
}

InitialBoundaryValueProblemPDE::~InitialBoundaryValueProblemPDE() {}

//void InitialBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension)
//{
//    _timeDimension = dimension;
//}

//const Dimension & InitialBoundaryValueProblemPDE::timeDimension() const
//{
//    return _timeDimension;
//}

//void InitialBoundaryValueProblemPDE::setSpaceDimensionX(const Dimension &dimension)
//{
//    _spaceDimensionX = dimension;
//}

//const Dimension& InitialBoundaryValueProblemPDE::spaceDimensionX() const
//{
//    return _spaceDimensionX;
//}

//void InitialBoundaryValueProblemPDE::setSpaceDimensionY(const Dimension &dimension)
//{
//    _spaceDimensionY = dimension;
//}

//const Dimension& InitialBoundaryValueProblemPDE::spaceDimensionY() const
//{
//    return _spaceDimensionY;
//}

//void InitialBoundaryValueProblemPDE::setSpaceDimensionZ(const Dimension &dimension)
//{
//    _spaceDimensionZ = dimension;
//}

//const Dimension& InitialBoundaryValueProblemPDE::spaceDimensionZ() const
//{
//    return _spaceDimensionZ;
//}

//void InitialBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY)
//{
//    this->_spaceDimensionX = dimensionX;
//    this->_spaceDimensionY = dimensionY;
//}

//void InitialBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension& dimensionZ)
//{
//    this->_spaceDimensionX = dimensionX;
//    this->_spaceDimensionY = dimensionY;
//    this->_spaceDimensionZ = dimensionZ;
//}

/*******************************************************************************************************************************************************/

FinalBoundaryValueProblemPDE::FinalBoundaryValueProblemPDE() {}

FinalBoundaryValueProblemPDE::FinalBoundaryValueProblemPDE(const FinalBoundaryValueProblemPDE &other) :
    FinalValueProblemPDE(other), BoundaryValueProblemPDE(other)/*,
    _timeDimension(other._timeDimension), _spaceDimensionX(other._spaceDimensionX),
    _spaceDimensionY(other._spaceDimensionY), _spaceDimensionZ(other._spaceDimensionZ)*/ {}

FinalBoundaryValueProblemPDE & FinalBoundaryValueProblemPDE::operator =(const FinalBoundaryValueProblemPDE &other)
{
    C_UNUSED(other);
//    this->_timeDimension = other._timeDimension;
//    this->_spaceDimensionX = other._spaceDimensionX;
//    this->_spaceDimensionY = other._spaceDimensionY;
//    this->_spaceDimensionZ = other._spaceDimensionZ;
    return *this;
}

FinalBoundaryValueProblemPDE::~FinalBoundaryValueProblemPDE() {}

//auto FinalBoundaryValueProblemPDE::setTimeDimension(const Dimension &dimension) -> void
//{
//    _timeDimension = dimension;
//}

//auto FinalBoundaryValueProblemPDE::timeDimension() const -> const Dimension &
//{
//    return _timeDimension;
//}

//auto FinalBoundaryValueProblemPDE::setSpaceDimensionX(const Dimension &dimension) -> void
//{
//    _spaceDimensionX = dimension;
//}

//auto FinalBoundaryValueProblemPDE::spaceDimensionX() const -> const Dimension&
//{
//    return _spaceDimensionX;
//}

//auto FinalBoundaryValueProblemPDE::setSpaceDimensionY(const Dimension &dimension) -> void
//{
//    _spaceDimensionY = dimension;
//}

//auto FinalBoundaryValueProblemPDE::spaceDimensionY() const -> const Dimension&
//{
//    return _spaceDimensionY;
//}

//auto FinalBoundaryValueProblemPDE::setSpaceDimensionZ(const Dimension &dimension) -> void
//{
//    _spaceDimensionZ = dimension;
//}

//auto FinalBoundaryValueProblemPDE::spaceDimensionZ() const -> const Dimension&
//{
//    return _spaceDimensionZ;
//}

//auto FinalBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY) -> void
//{
//    this->_spaceDimensionX = dimensionX;
//    this->_spaceDimensionY = dimensionY;
//}

//auto FinalBoundaryValueProblemPDE::setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension& dimensionZ) -> void
//{
//    this->_spaceDimensionX = dimensionX;
//    this->_spaceDimensionY = dimensionY;
//    this->_spaceDimensionZ = dimensionZ;
//}
