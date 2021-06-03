#ifndef INITIAL_BOUNDARY_VALUE_PROBLEM_H
#define INITIAL_BOUNDARY_VALUE_PROBLEM_H

#include "ivp.h"
#include "bvp.h"

class MINIMUMSHARED_EXPORT GridDimensionPDE
{
protected:
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(GridDimensionPDE)

    virtual Dimension timeDimension() const = 0;
    virtual Dimension spaceDimensionX() const = 0;
    virtual Dimension spaceDimensionY() const = 0;
    virtual Dimension spaceDimensionZ() const = 0;
};

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemODE : public InitialValueProblemODE, public BoundaryValueProblemODE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialBoundaryValueProblemODE)

public:
    virtual Dimension dimension() const = 0;
};

class MINIMUMSHARED_EXPORT FinalBoundaryValueProblemODE : public FinalValueProblemODE, public BoundaryValueProblemODE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalBoundaryValueProblemODE)

public:
    virtual Dimension dimension() const = 0;
};

class MINIMUMSHARED_EXPORT InitialBoundaryValueProblemPDE : public InitialValueProblemPDE, public BoundaryValueProblemPDE, public GridDimensionPDE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(InitialBoundaryValueProblemPDE)
};

class MINIMUMSHARED_EXPORT FinalBoundaryValueProblemPDE : public FinalValueProblemPDE, public BoundaryValueProblemPDE, public GridDimensionPDE
{
    PUBLIC_CONSTRUCTORS_VIRTUAL_DESTRUCTOR(FinalBoundaryValueProblemPDE)
};

#endif // INITIAL_BOUNDARY_VALUE_PROBLEM_H
