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
class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : public InitialValueProblemPDE, public BoundaryValueProblemPDE
{
public:
    InitialBoundaryValueProblemPDE();
    InitialBoundaryValueProblemPDE(const InitialBoundaryValueProblemPDE &);
    InitialBoundaryValueProblemPDE & operator = (const InitialBoundaryValueProblemPDE &);
    virtual ~InitialBoundaryValueProblemPDE();

    virtual const Dimension& timeDimension() const;
    virtual const Dimension& spaceDimensionX() const;
    virtual const Dimension& spaceDimensionY() const;
    virtual const Dimension& spaceDimensionZ() const;

    virtual void setTimeDimension(const Dimension &dimension);
    virtual void setSpaceDimensionX(const Dimension &dimensionX);
    virtual void setSpaceDimensionY(const Dimension &dimensionY);
    virtual void setSpaceDimensionZ(const Dimension &dimensionZ);

    virtual void setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY);
    virtual void setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension& dimensionZ);

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};

/**
 * @brief The InitialBoundaryValueProblemPDE class
 * @see ParabolicIBVP
 */
class MINIMUMSHARED_EXPORT FinalBoundaryValueProblemPDE : public FinalValueProblemPDE, public BoundaryValueProblemPDE
{
public:
    FinalBoundaryValueProblemPDE();
    FinalBoundaryValueProblemPDE(const FinalBoundaryValueProblemPDE &);
    FinalBoundaryValueProblemPDE & operator = (const FinalBoundaryValueProblemPDE &);
    virtual ~FinalBoundaryValueProblemPDE();

    virtual const Dimension& timeDimension() const;
    virtual const Dimension& spaceDimensionX() const;
    virtual const Dimension& spaceDimensionY() const;
    virtual const Dimension& spaceDimensionZ() const;

    virtual void setTimeDimension(const Dimension &dimension);
    virtual void setSpaceDimensionX(const Dimension &dimension);
    virtual void setSpaceDimensionY(const Dimension &dimension);
    virtual void setSpaceDimensionZ(const Dimension &dimension);

    virtual void setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY);
    virtual void setSpaceDimensions(const Dimension& dimensionX, const Dimension& dimensionY, const Dimension& dimensionZ);

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
