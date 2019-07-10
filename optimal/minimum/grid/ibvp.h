#ifndef INITIAL_BOUNDARY_VALUE_PROBLEM_H
#define INITIAL_BOUNDARY_VALUE_PROBLEM_H

#include "grid.h"
#include "ivp.h"
#include "bvp.h"
#include <math.h>
#include "../cmethods.h"
#include "../linearequation.h"
#include "../matrix3d.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

/**
 * @brief The InitialBoundaryValueProblemPDE class
 * @see ParabolicIBVP
 */
class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE :
        public InitialValueProblemPDE,
        public BoundaryValueProblemPDE
{
public:
    InitialBoundaryValueProblemPDE();
    InitialBoundaryValueProblemPDE(const InitialBoundaryValueProblemPDE& ibvp);
    InitialBoundaryValueProblemPDE& operator =(const InitialBoundaryValueProblemPDE&);
    virtual ~InitialBoundaryValueProblemPDE();

    virtual auto setTimeDimension(const Dimension &dimension) -> void;
    virtual auto timeDimension() const -> const Dimension&;
    virtual auto setSpaceDimensionX(const Dimension &dimension) -> void;
    virtual auto spaceDimensionX() const -> const Dimension&;
    virtual auto setSpaceDimensionY(const Dimension &dimension) -> void;
    virtual auto spaceDimensionY() const -> const Dimension&;
    virtual auto setSpaceDimensionZ(const Dimension &dimension) -> void;
    virtual auto spaceDimensionZ() const -> const Dimension&;

    //virtual auto addSpaceDimension(const Dimension &dimension) -> void;
    //virtual auto spaceDimension(Dimension::SpaceDimension dim) const -> const Dimension&;
    //virtual auto dimSize() const -> unsigned int;

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
    //std::vector<Dimension> mspaceDimension;
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
