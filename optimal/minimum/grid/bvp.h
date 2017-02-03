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
        Right = 1,
        Unused = 2
    };

    BoundaryCondition condition;
};

struct MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
    virtual double boundary(BoundaryType bound) const = 0;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
};

#endif // BOUNDARYVALUEPROBLEM_H
