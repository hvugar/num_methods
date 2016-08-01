#ifndef PROJECTION_H
#define PROJECTION_H

#include "global.h"
#include <vector2d.h>

class MINIMUMSHARED_EXPORT IProjection
{
public:
    virtual void project(DoubleVector &x, int index) = 0;
};

#endif // PROJECTION_H
