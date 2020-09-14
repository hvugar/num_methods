#ifndef REGION_ELIMINATION_METHOD_H
#define REGION_ELIMINATION_METHOD_H

#include "global.h"
#include "vector2d.h"

/**
 * @brief The Abstract Gradient Method class
 */
class MINIMUMSHARED_EXPORT RegionEliminationMethod
{
public:
    RegionEliminationMethod();
    virtual ~RegionEliminationMethod();

    virtual void calculate(DoubleVector &x) = 0;
};

#endif // REGION_ELIMINATION_METHOD_H
