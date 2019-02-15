#include "integral1.h"

NumericalIntegral::NumericalIntegral(const Dimension &dimension) : _dimension(dimension)
{
}

double NumericalIntegral::rectangleRule(RectangleRuleType type) const
{
    const int min = dimension().min();
    const int max = dimension().max();
    const double h = dimension().step();

    double sum = 0.0;

    if (type == Left)
    {
        for (int n=min; n<max; n++)
        {
            sum += f(n*h, n);
        }
    }

    if (type == Right)
    {
        for (int n=min+1; n<=max; n++)
        {
            sum += f(n*h, n);
        }
    }

    if (type == Center)
    {
        for (int n=min; n<max; n++)
        {
            sum += f( 0.5*(n*h+(n+1)*h), n);
        }
    }
    sum *= h;
    return sum;
}

double NumericalIntegral::trapezoidalRule() const
{
    const int min = dimension().min();
    const int max = dimension().max();
    const double h = dimension().step();

    double sum = 0.0;
    for (int n=min; n<max; n++)
    {
        sum += ( f(n*h, n) + f((n+1)*h, n+1) );
    }
    sum *= 0.5*h;
    return sum;
}

double NumericalIntegral::SimpsonsRule() const
{
    return 0.0;
}

const Dimension& NumericalIntegral::dimension() const
{
    return _dimension;
}

void NumericalIntegral::setDimension(const Dimension& dimension)
{
    _dimension = dimension;
}
