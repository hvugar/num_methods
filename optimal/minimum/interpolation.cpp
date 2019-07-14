#include "interpolation.h"

Interpolation::Interpolation()
{}

LagrangeInterpolation2D::~LagrangeInterpolation2D()
{}

double LagrangeInterpolation2D::value(double, double, const DoubleMatrix &) const
{
    return 0.0;
}
