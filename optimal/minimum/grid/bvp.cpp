#include "bvp.h"

const UniformGrid& BoundaryValueProblem::spaceGridX() const
{
    return spacegridX;
}

const UniformGrid& BoundaryValueProblem::spaceGridY() const
{
    return spacegridY;
}

const UniformGrid& BoundaryValueProblem::spaceGridZ() const
{
    return spacegridZ;
}

void BoundaryValueProblem::setSpaceGridX(const UniformGrid& spaceGridX)
{
    spacegridX = spaceGridX;
}

void BoundaryValueProblem::setSpaceGridY(const UniformGrid& spaceGridY)
{
    spacegridY = spaceGridY;
}

void BoundaryValueProblem::setSpaceGridZ(const UniformGrid& spaceGridZ)
{
    spacegridZ = spaceGridZ;
}

