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
    virtual ~InitialBoundaryValueProblemPDE();

    auto setTimeDimension(const Dimension &dimension) -> void;
    auto timeDimension() const -> const Dimension&;

    auto addSpaceDimension(const Dimension &dimension) -> void;
    auto spaceDimension(Dimension::SpaceDimension dim) const -> const Dimension&;

    auto dimSize() const -> unsigned int;

protected:
    Dimension mtimeDimension;
    std::vector<Dimension> mspaceDimension;
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
