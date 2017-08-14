#include "numintegralexp1.h"

void NumIntegralExp1::Main(int agrc, char *argv[])
{
    ODEGrid grid(Dimension(0.01, 100, 0));
    NumIntegralExp1 exp1(grid);
    double rr = exp1.rectangleRule(Left);
    printf("%f\n", rr);
}

NumIntegralExp1::NumIntegralExp1(const ODEGrid &grid) : NumericalIntegral(grid)
{

}

double NumIntegralExp1::f(double x, int n UNUSED_PARAM) const
{
    return 2.0*x;
}
