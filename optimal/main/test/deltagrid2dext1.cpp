#include "deltagrid2dext1.h"

void DeltaGrid2DExt1::Main(int argc, char **argv)
{
    DeltaGrid2DExt1 ex1;

    unsigned int N = 100;
    unsigned int M = 100;
    DoubleMatrix mx(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            mx[m][n] = ex1.fx(0.01*n, 0.01*m);
        }
    }

    DeltaGrid2D dg;
    dg.initGrid(100, 0.01, 100, 0.01);
    dg.distributeGauss(SpacePoint(0.4224, 0.6538));

    double dx, dy;
    double z = dg.consentrateInPoint(mx, dx, dy);

    printf("%f %f %f %f %f\n", z, dx, dy, dg.p().x, dg.p().y);
    printf("%f %f %f\n", ex1.fx(dg.p().x, dg.p().y),
                         ex1.dx(dg.p().x, dg.p().y),
                         ex1.dy(dg.p().x, dg.p().y));
}

double DeltaGrid2DExt1::fx(double x, double y) const
{
    return sin(8.0*M_PI*x)*sin(8.0*M_PI*y);
}

double DeltaGrid2DExt1::dx(double x, double y) const
{
    return 8.0*M_PI*cos(8.0*M_PI*x)*sin(8.0*M_PI*y);
}

double DeltaGrid2DExt1::dy(double x, double y) const
{
    return 8.0*M_PI*sin(8.0*M_PI*x)*cos(8.0*M_PI*y);
}
