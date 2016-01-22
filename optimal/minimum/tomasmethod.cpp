#include "tomasmethod.h"
#include <cmethods.h>

void TomasAlgorithm(const DoubleVector &a, const DoubleVector &b, const DoubleVector &c, const DoubleVector &d, DoubleVector &x)
{
    tomasAlgorithm(a.data(), b.data(), c.data(), d.data(), x.data(), x.size());
}
