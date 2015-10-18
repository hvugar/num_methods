#ifndef PROJECTION_H
#define PROJECTION_H

#include "global.h"
#include "doublevector.h"

class MINIMUMSHARED_EXPORT Projection
{
public:
    virtual void project(DoubleVector &x, int index);

    double a;
    double b;
};

#endif // PROJECTION_H
