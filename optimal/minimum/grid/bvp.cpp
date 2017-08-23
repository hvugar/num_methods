#include "bvp.h"

const UniformODEGrid& BoundaryValueProblem::spaceGridX() const
{
    return spacegridX;
}

const UniformODEGrid& BoundaryValueProblem::spaceGridY() const
{
    return spacegridY;
}

const UniformODEGrid& BoundaryValueProblem::spaceGridZ() const
{
    return spacegridZ;
}

void BoundaryValueProblem::setSpaceGridX(const UniformODEGrid& spaceGridX)
{
    spacegridX = spaceGridX;
}

void BoundaryValueProblem::setSpaceGridY(const UniformODEGrid& spaceGridY)
{
    spacegridY = spaceGridY;
}

void BoundaryValueProblem::setSpaceGridZ(const UniformODEGrid& spaceGridZ)
{
    spacegridZ = spaceGridZ;
}

