#ifndef PROJECTION_H
#define PROJECTION_H

#include "doublevector.h"

class Projection
{
public:
    virtual void project(DoubleVector &x, int index);

    double a;
    double b;
};

#endif // PROJECTION_H
