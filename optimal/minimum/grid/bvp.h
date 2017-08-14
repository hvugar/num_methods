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
    const UniformGrid& spaceGridX() const;
    const UniformGrid& spaceGridY() const;
    const UniformGrid& spaceGridZ() const;

    void setSpaceGridX(const UniformGrid&);
    void setSpaceGridY(const UniformGrid&);
    void setSpaceGridZ(const UniformGrid&);

private:
    UniformGrid spacegridX;
    UniformGrid spacegridY;
    UniformGrid spacegridZ;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemODE : protected BoundaryValueProblem
{
protected:
    virtual double boundary(BoundaryType bound) const = 0;
};

class MINIMUMSHARED_EXPORT BoundaryValueProblemPDE : protected BoundaryValueProblem
{
};

#endif // BOUNDARYVALUEPROBLEM_H
