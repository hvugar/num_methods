#ifndef BOUNDARYVALUEPROBLEM_H
#define BOUNDARYVALUEPROBLEM_H

#include "../vector2d.h"
#include "../matrix2d.h"
#include "../matrix2d.h"
#include "../cmethods.h"
#include "../printer.h"
#include <math.h>
#include "grid.h"

struct MINIMUMSHARED_EXPORT InitialValueProblem
{
    Dimension timeDimension;
};

struct MINIMUMSHARED_EXPORT BoundaryValueProblem
{
    enum BoundaryCondition
    {
        Dirichlet = 1,
        Neumann = 2,
        Mixed = 3
    };

    enum BoundaryType
    {
        Left = 0,
        Right = 1
    };

    BoundaryCondition condition;
};

struct MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
    virtual double boundary(BoundaryType bound) const = 0;
};

struct MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
    virtual double boundary(unsigned int m, BoundaryType bound) const = 0;

    SpaceDimension spaceDimension;
};

#endif // BOUNDARYVALUEPROBLEM_H
