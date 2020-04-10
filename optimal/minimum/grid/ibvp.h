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

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemODE : public InitialValueProblemODE, public BoundaryValueProblemODE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialBoundaryValueProblemODE);

public:
    virtual Dimension dimension() const = 0;
};

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : public InitialValueProblemPDE, public BoundaryValueProblemPDE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialBoundaryValueProblemPDE);

public:
    virtual Dimension timeDimension() const = 0;
    virtual Dimension spaceDimensionX() const = 0;
    virtual Dimension spaceDimensionY() const = 0;
    virtual Dimension spaceDimensionZ() const = 0;
};

class MINIMUMSHARED_EXPORT FinalBoundaryValueProblemPDE : public FinalValueProblemPDE, public BoundaryValueProblemPDE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalBoundaryValueProblemPDE);

public:
    virtual Dimension timeDimension() const = 0;
    virtual Dimension spaceDimensionX() const = 0;
    virtual Dimension spaceDimensionY() const = 0;
    virtual Dimension spaceDimensionZ() const = 0;
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
