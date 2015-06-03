#include "gradient.h"

Gradient::~Gradient() {}

double Gradient::epsilon() const
{
    return mepsilon;
}

void Gradient::setEpsilon(double epsilon)
{
    mepsilon = epsilon;
}

void Gradient::setX(const std::vector<double> &x)
{
    mx = x;
}

const std::vector<double>& Gradient::x() const
{
    return mx;
}
