#include "integral1.h"

NumericalIntegral::NumericalIntegral(const ODEGrid &grid) : mgrid(grid)
{
}

double NumericalIntegral::rectangleRule(RectangleRuleType type) const
{
    Dimension dim = grid().dimension();
    unsigned int minN = dim.minN();
    unsigned int maxN = dim.maxN();
    double h = dim.step();

    double sum = 0.0;
    if (type == Left)
    {
        for (unsigned int n=minN; n<maxN; n++)
        {
            sum += f(n*h, n);
        }
    }
    if (type == Right)
    {
        for (unsigned int n=minN+1; n<=maxN; n++)
        {
            sum += f(n*h, n);
        }
    }
    if (type == Center)
    {
        for (unsigned int n=minN; n<maxN; n++)
        {
            sum += f( 0.5*(n*h+(n+1)*h), n);
        }
    }
    sum *= h;
    return sum;
}

double NumericalIntegral::trapezoidalRule() const
{
    Dimension dim = grid().dimension();
    unsigned int minN = dim.minN();
    unsigned int maxN = dim.maxN();
    double h = dim.step();

    double sum = 0.0;
    for (unsigned int n=minN; n<maxN; n++)
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

const ODEGrid& NumericalIntegral::grid() const
{
    return mgrid;
}

void NumericalIntegral::setGrid(const ODEGrid& grid)
{
    mgrid = grid;
}
