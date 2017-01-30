#include "parabolicibvp1.h"

double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
{
    TimeDimension time = grid().timeDimension();
    SpaceDimension space = grid().spaceDimensions(SpaceDimension::Dim1);

    double x = n * space.hx();
    double t = m * time.ht();
    return x*x + t;
}

double ParabolicIBVP1::initial(unsigned int n) const
{
    return U(n, 0);
}

double ParabolicIBVP1::boundary(unsigned int m, BoundaryType boundary) const
{
    SpaceDimension space = grid().spaceDimensions(SpaceDimension::Dim1);
    unsigned int N = space.N2();

    if (boundary == Left) return U(0,m);
    if (boundary == Right) return U(N,m);
    return 0.0;
}

double ParabolicIBVP1::f(unsigned int n, unsigned int m) const
{
    return 1.0 - 2.0*a(n,m);
}

double ParabolicIBVP1::a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const
{
    return 1.0;
}
