#include "bvp.h"

const UniformODEGrid& BoundaryValueProblem::spaceGridX() const
{
    return mSpacegridX;
}

const UniformODEGrid& BoundaryValueProblem::spaceGridY() const
{
    return mSpacegridY;
}

const UniformODEGrid& BoundaryValueProblem::spaceGridZ() const
{
    return mSpacegridZ;
}

void BoundaryValueProblem::setSpaceGridX(const UniformODEGrid& spaceGridX)
{
    mSpacegridX = spaceGridX;
}

void BoundaryValueProblem::setSpaceGridY(const UniformODEGrid& spaceGridY)
{
    mSpacegridY = spaceGridY;
}

void BoundaryValueProblem::setSpaceGridZ(const UniformODEGrid& spaceGridZ)
{
    mSpacegridZ = spaceGridZ;
}

