#include "deltagrid2dext1.h"

void DeltaGrid2DExt1::Main(int argc, char **argv)
{
    DeltaGrid2DExt1 ex1;

    unsigned int N = 100; double hx = 0.01;
    unsigned int M = 100; double hy = 0.01;
    DoubleMatrix mx(M+1, N+1);
    for (unsigned int m=0; m<=M; m++)
    {
        for (unsigned int n=0; n<=N; n++)
        {
            mx[m][n] = ex1.fx(hx*n, hy*m);
            //mx[m][n] = 1.0;
        }
    }

    DeltaGrid2D dg;
    dg.initGrid(N, hx, M, hy);
    dg.distributeGauss(SpacePoint(0.421, 0.651), 1, 1);
    //dg.distributeGauss(SpacePoint(0.4224, 0.6538), 1, 1);

    double dx, dy;
    double z = dg.consentrateInPoint(mx, dx, dy);
    //double z = dg.consentrateInPoint(mx, 4);

    printf("%10.6f %10.6f %10.6f\n", ex1.fx(dg.p().x, dg.p().y), ex1.dx(dg.p().x, dg.p().y), ex1.dy(dg.p().x, dg.p().y));
    //printf("%10.6f %10.6f %10.6f\n", ex1.fx(dg.rx()*dg.hx(), dg.ry()*dg.hy()), dg.rx()*1.0, dg.ry()*1.0);
    puts("---");
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", dg.consentrateInPoint(mx, 1), dx, dy, dg.p().x, dg.p().y);
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", dg.consentrateInPoint(mx, 2), dx, dy, dg.p().x, dg.p().y);
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", dg.consentrateInPoint(mx, 0), dx, dy, dg.p().x, dg.p().y);
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", dg.consentrateInPoint(mx, 3), dx, dy, dg.p().x, dg.p().y);
    printf("%10.6f %10.6f %10.6f %10.6f %10.6f\n", dg.consentrateInPoint(mx, 4), dx, dy, dg.p().x, dg.p().y);
}

double DeltaGrid2DExt1::fx(double x, double y) const
{
    //return sin(8.0*M_PI*x)*sin(8.0*M_PI*y);
    return sin(50.0*x*x)*cos(5.0*x)+sin(20.0*y*y)*sin(5.0*y)-0.5;
}

double DeltaGrid2DExt1::dx(double x, double y) const
{
    //return 8.0*M_PI*cos(8.0*M_PI*x)*sin(8.0*M_PI*y);
    return 100.0*x*cos(50.0*x*x)*cos(5.0*x) - 5.0*sin(50.0*x*x)*sin(5.0*x);
}

double DeltaGrid2DExt1::dy(double x, double y) const
{
    //return 8.0*M_PI*sin(8.0*M_PI*x)*cos(8.0*M_PI*y);
    return 40.0*y*cos(20.0*y*y)*sin(5.0*y) + 5.0*sin(20.0*y*y)*cos(5.0*y);
}
