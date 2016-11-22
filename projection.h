#ifndef PROJECTION_H
#define PROJECTION_H

#include "doublevector.h"

class IProjection
{
public:
    virtual void project(DoubleVector &x, int index) = 0;
};

#endif // PROJECTION_H
