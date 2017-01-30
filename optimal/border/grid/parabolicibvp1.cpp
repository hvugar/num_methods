#include "parabolicibvp1.h"

void ParabolicIBVP1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    GridPDE grid;
    grid.setTimeDimension(TimeDimension(0.001, 0.0, 1.0, 0, 1000));
    grid.addSpaceDimension(SpaceDimension(0.001, 0.0, 1.0, 0, 1000));

    ParabolicIBVP1 p(grid);
    p.setGrid(grid);
    DoubleMatrix u0;
    p.gridMethod(u0);
    IPrinter::printMatrix(14,10,u0);
    IPrinter::printSeperatorLine();
    DoubleMatrix u1;
    p.calculateN2L2RD(u1);
    IPrinter::printMatrix(14,10,u1);
    IPrinter::printSeperatorLine();
    DoubleMatrix u2;
    p.calculateN4L2RD(u2);
    IPrinter::printMatrix(14,10,u2);
}

ParabolicIBVP1::ParabolicIBVP1(const GridPDE &grid) : ParabolicIBVP()
{
    setGrid(grid);
}

double ParabolicIBVP1::initial(unsigned int n) const
{
    return U(n,grid().timeDimension().M1());
}

double ParabolicIBVP1::boundary(unsigned int m, BoundaryType boundary) const
{
    SpaceDimension dim1 =  grid().spaceDimensions(SpaceDimension::Dim1);
    if (boundary == Left)  return U(0,m);
    if (boundary == Right) return U(dim1.N2(),m);
    return 0.0;
}

double ParabolicIBVP1::f(unsigned int n, unsigned int m) const
{
    double t = m*grid().timeDimension().ht();
    //double x = n*grid().spaceDimensions(SpaceDimension::Dim1).hx();

#ifdef SAMPLE_1
    return 1.0 - 2.0*a(n,m);
#endif
#ifdef SAMPLE_2
    return 2.0*t - 2.0*a(n,m);
#endif
    return 0.0;
}

double ParabolicIBVP1::a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*grid().timeDimension().ht();
    double x = n*grid().spaceDimensions(SpaceDimension::Dim1).hx();

#ifdef SAMPLE_1
    return x*x + t;
#endif
#ifdef SAMPLE_2
    return x*x + t*t;
#endif
    return 0.0;
}
