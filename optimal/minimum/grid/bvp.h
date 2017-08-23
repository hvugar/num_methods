#ifndef BOUNDARYVALUEPROBLEM_H
#define BOUNDARYVALUEPROBLEM_H

#include "../vector2d.h"
#include "../matrix2d.h"
#include "../matrix3d.h"
#include "../cmethods.h"
#include "../printer.h"
#include <math.h>
#include "grid.h"
#include "uniformgrid.h"

class MINIMUMSHARED_EXPORT BoundaryValueProblem
{
public:
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

public:
    const UniformODEGrid& spaceGridX() const;
    const UniformODEGrid& spaceGridY() const;
    const UniformODEGrid& spaceGridZ() const;

    void setSpaceGridX(const UniformODEGrid&);
    void setSpaceGridY(const UniformODEGrid&);
    void setSpaceGridZ(const UniformODEGrid&);

private:
    UniformODEGrid spacegridX;
    UniformODEGrid spacegridY;
    UniformODEGrid spacegridZ;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
protected:
    virtual double boundary(BoundaryType bound) const = 0;
};

/**
 * @brief The BoundaryValueProblemPDE class
 * @see
 */
class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
public:
    double alphaL;
    double alphaR;
    double bettaL;
    double bettaR;

    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
};

#endif // BOUNDARYVALUEPROBLEM_H
