#include "projection.h"

void Projection::project(DoubleVector &x)
{
}

void Projection::project(double &x)
{
    if (x < a) x = a;
    else
    if (x > b) x = b;
}
