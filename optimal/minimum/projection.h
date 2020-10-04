#ifndef PROJECTION_H
#define PROJECTION_H

#include "global.h"
#include "vector2d.h"

class MINIMUMSHARED_EXPORT IProjection
{
public:
    virtual ~IProjection();

    virtual void project(DoubleVector &x, size_t index) = 0;
    virtual void project(DoubleVector &) const {}
};

#endif // PROJECTION_H
