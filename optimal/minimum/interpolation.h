#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "vector2d.h"
#include "matrix2d.h"

class Interpolation
{
public:
    Interpolation();
};

class LagrangeInterpolation1D
{};

class LagrangeInterpolation2D : public Interpolation
{
public:
    virtual ~LagrangeInterpolation2D();

    virtual double value(double x, double y, const DoubleMatrix &m) const;

private:
    unsigned int _Nx;
    unsigned int _Ny;
    double _hx;
    double _hy;
    double _minX;
    double _maxX;
    double _minY;
    double _maxY;
};

#endif // INTERPOLATION_H
