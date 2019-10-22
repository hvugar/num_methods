#ifndef INITIAL_BOUNDARY_VALUE_PROBLEM_H
#define INITIAL_BOUNDARY_VALUE_PROBLEM_H

#include "ivp.h"
#include "bvp.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

/**
 * @brief The InitialBoundaryValueProblemPDE class
 * @class InitialBoundaryValueProblemPDE
 * @see InitialValueProblemPDE
 * @see BoundaryValueProblemPDE
 * @copyright
 */
class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : public InitialValueProblemPDE, public BoundaryValueProblemPDE
{
public:
    InitialBoundaryValueProblemPDE();
    InitialBoundaryValueProblemPDE(const InitialBoundaryValueProblemPDE &);
    InitialBoundaryValueProblemPDE & operator= (const InitialBoundaryValueProblemPDE &);
    virtual ~InitialBoundaryValueProblemPDE();

    virtual const Dimension& timeDimension() const = 0;
    virtual const Dimension& spaceDimensionX() const = 0;
    virtual const Dimension& spaceDimensionY() const = 0;
    virtual const Dimension& spaceDimensionZ() const = 0;
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

    virtual const Dimension& timeDimension() const = 0;
    virtual const Dimension& spaceDimensionX() const = 0;
    virtual const Dimension& spaceDimensionY() const = 0;
    virtual const Dimension& spaceDimensionZ() const = 0;
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
