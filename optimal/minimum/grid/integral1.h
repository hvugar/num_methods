#ifndef INTEGRAL1_H
#define INTEGRAL1_H

#include "grid.h"

class MINIMUMSHARED_EXPORT NumericalIntegral
{
public:
    enum RectangleRuleType
    {
        Left,
        Center,
        Right
    };

    NumericalIntegral(const UniformODEGrid &grid);

    const UniformODEGrid& grid() const;
    void setGrid(const UniformODEGrid& grid);

    virtual double f(double x, int n) const = 0;

    double rectangleRule(RectangleRuleType type = Left) const;
    double trapezoidalRule() const;
    double SimpsonsRule() const;


private:
    UniformODEGrid mgrid;
};

#endif // INTEGRAL1_H
