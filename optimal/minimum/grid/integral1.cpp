#include "integral1.h"

NumericalIntegral::NumericalIntegral(const UniformODEGrid &grid) : mgrid(grid)
{
}

double NumericalIntegral::rectangleRule(RectangleRuleType type) const
{
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const double h = grid().dimension().step();

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
    const int min = grid().dimension().min();
    const int max = grid().dimension().max();
    const double h = grid().dimension().step();

    double sum = 0.0;
    for (unsigned int n=min; n<max; n++)
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

const UniformODEGrid& NumericalIntegral::grid() const
{
    return mgrid;
}

void NumericalIntegral::setGrid(const UniformODEGrid& grid)
{
    mgrid = grid;
}
