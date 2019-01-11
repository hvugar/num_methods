#ifndef BOUNDARYVALUEPROBLEM_H
#define BOUNDARYVALUEPROBLEM_H

#include "../vector2d.h"
#include "../matrix2d.h"
#include "../matrix3d.h"
#include "../cmethods.h"
#include "../printer.h"
#include "../linearequation.h"
#include <math.h>
#include "grid.h"

class MINIMUMSHARED_EXPORT BoundaryValueProblem
{
public:
    enum BoundaryCondition
    {
        Dirichlet = 1,
        Neumann = 2,
        Mixed = 3
    };

    BoundaryCondition condition;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
public:
    virtual ~BoundaryValueProblemODE();

protected:
    virtual double boundary() const = 0;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
public:
    virtual ~BoundaryValueProblemPDE();

protected:
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
};

#endif // BOUNDARYVALUEPROBLEM_H
