#include "projection.h"

void Projection::project(DoubleVector &x, int index)
{
    if (x[index] < a) x[index] = a;
    else
    if (x[index] > b) x[index] = b;
}
