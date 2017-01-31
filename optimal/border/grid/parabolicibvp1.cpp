#include "parabolicibvp1.h"
#include <time.h>

void ParabolicIBVP1::Main(int agrc UNUSED_PARAM, char *argv[] UNUSED_PARAM)
{
    GridPDE grid;
    grid.setTimeDimension(TimeDimension(0.0001, 0.0, 1.0, 0, 10000));
    grid.addSpaceDimension(SpaceDimension(0.0001, 0.0, 1.0, 0, 10000));

    ParabolicIBVP1 p(grid);
    clock_t t;
    DoubleMatrix u0;
    t = clock();
    p.gridMethod(u0);
    t = clock() - t;
    IPrinter::printMatrix(14,10,u0);
    printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    IPrinter::printSeperatorLine();

//    DoubleMatrix u1;
//    p.calculateN2L2RD(u1);
//    IPrinter::printMatrix(14,10,u1);
//    IPrinter::printSeperatorLine();

//    DoubleMatrix u2;
//    p.calculateN4L2RD(u2);
//    IPrinter::printMatrix(14,10,u2);
//    IPrinter::printSeperatorLine();

    DoubleMatrix u01;
    t = clock();
    p.gridMethod1(u01);
    t = clock() - t;
    IPrinter::printMatrix(14,10,u01);
    printf ("It took me %ld clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
    IPrinter::printSeperatorLine();
}

ParabolicIBVP1::ParabolicIBVP1(const GridPDE &grid) : ParabolicIBVP()
{
    setGrid(grid);
}

double ParabolicIBVP1::initial(unsigned int n) const
{
    return U(n,grid().timeDimension().minN());
}

double ParabolicIBVP1::initial(const SpaceNode& n) const
{
    return U(n.i,grid().timeDimension().minN());
}

double ParabolicIBVP1::boundary(unsigned int m, BoundaryType boundary) const
{
    SpaceDimension dim1 =  grid().spaceDimensions(SpaceDimension::Dim1);
    if (boundary == Left)  return U(dim1.minN(),m);
    if (boundary == Right) return U(dim1.maxN(),m);
    return 0.0;
}

double ParabolicIBVP1::boundary(const SpaceNode &sn, const TimeNode &tn) const
{
    return U(sn.i,tn.i);
}

double f1(double x, double t, double a)
{
#ifdef SAMPLE_1
    return 1.0 - 2.0*a;
#endif
#ifdef SAMPLE_2
    return 2.0*t - 2.0*a;
#endif
#ifdef SAMPLE_3
    return sin(x)*(1.0+a*t);
#endif
#ifdef SAMPLE_4
    return (sin(10.0*x) + 2.0*x*exp(2.0*x*t) - a*(4.0*t*t*exp(2.0*x*t) - 100.0*t*sin(10.0*x)));
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*(1.0 - a*k*k*t);
#endif
}

double ParabolicIBVP1::f(unsigned int n, unsigned int m) const
{
    double t = m*grid().timeDimension().step();
    double x = n*grid().spaceDimensions(SpaceDimension::Dim1).step();
    return f1(x,t,a(n,m));
}

double ParabolicIBVP1::f(const SpaceNode &sn, const TimeNode &tn) const
{
    double x = sn.x;
    double t = tn.t;
    return f1(x,t,a(sn,tn));
}

double ParabolicIBVP1::a(unsigned int n UNUSED_PARAM, unsigned int m UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP1::a(const SpaceNode &sn UNUSED_PARAM, const TimeNode &tn UNUSED_PARAM) const
{
    return 1.0;
}

double ParabolicIBVP1::U(unsigned int n, unsigned int m) const
{
    double t = m*grid().timeDimension().step();
    double x = n*grid().spaceDimensions(SpaceDimension::Dim1).step();

#ifdef SAMPLE_1
    return x*x + t;
#endif
#ifdef SAMPLE_2
    return x*x + t*t;
#endif
#ifdef SAMPLE_3
    return sin(x)*t;
#endif
#ifdef SAMPLE_4
    return t*sin(10.0*x) + exp(2.0*x*t);
#endif
#ifdef SAMPLE_5
    double k = 20.0;
    return exp(k*x)*t;
#endif
    return 0.0;
}
