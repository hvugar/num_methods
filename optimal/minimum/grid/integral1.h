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

    NumericalIntegral(const Dimension &dimension);

    const Dimension& dimension() const;
    void setDimension(const Dimension& dimension);

    virtual double f(double x, int n) const = 0;

    double rectangleRule(RectangleRuleType type = Left) const;
    double trapezoidalRule() const;
    double SimpsonsRule() const;


private:
    Dimension _dimension;
};

#endif // INTEGRAL1_H
