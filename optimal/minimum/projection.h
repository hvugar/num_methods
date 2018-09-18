#ifndef PROJECTION_H
#define PROJECTION_H

#include "global.h"
#include "vector2d.h"

class MINIMUMSHARED_EXPORT IProjection
{
public:
    virtual void project(DoubleVector &x, unsigned int index) = 0;
    virtual void project(DoubleVector &x) const {}
};

#endif // PROJECTION_H
