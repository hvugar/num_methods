#include "iproblem2pibvp2d.h"

double IProblem22DPIBVP::boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryType) const
{
    return NAN;
}

double IProblem22DPIBVP::f(const SpaceNodePDE&, const TimeNodePDE&) const
{
    return 0.0;
}

void IProblem22DPIBVP::setEquationParameters(double a, double lambda0, double lambda, double theta)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
    this->theta = theta;
}

void IProblem22DPIBVP::setSetting(P2Setting setting)
{
    this->setting = setting;
}

//-------------------------------------------------------------------------------------------------------//

ExtendedSpaceNode2D::ExtendedSpaceNode2D()
{
    x = y = 0.0;
    i = j = 0;
    wi = NULL;
    layerNumber = 0;
}

void ExtendedSpaceNode2D::setSpaceNode(const SpaceNodePDE &sn)
{
    this->i = sn.i;
    this->j = sn.j;
    this->x = sn.x;
    this->y = sn.y;
}

void ExtendedSpaceNode2D::extendWeights(const Dimension &dimX, const Dimension &dimY, unsigned int rows, unsigned int cols)
{
    this->rows = rows;
    this->cols = cols;
    wi = new WISpaceNodePDE*[rows];
    for (unsigned int j=0; j<rows; j++) wi[j] = new WISpaceNodePDE[cols];

    unsigned int Nx = dimX.sizeN();
    unsigned int Ny = dimY.sizeN();

    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int rx = (unsigned int)(floor(x*Nx));
    unsigned int ry = (unsigned int)(floor(y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = (1.0/(2.0*hy3));
    double hy36 = (1.0/(6.0*hy3));

    double dx = 0.0;
    double dy = 0.0;

    wi[1][1].i = rx + 0; wi[1][1].x = wi[1][1].i*hx;
    wi[1][1].j = ry + 0; wi[1][1].y = wi[1][1].j*hy;
    dx = fabs(wi[1][1].x-x);
    dy = fabs(wi[1][1].y-y);
    wi[1][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[2][1].i = rx + 0; wi[2][1].x = wi[2][1].i*hx;
    wi[2][1].j = ry + 1; wi[2][1].y = wi[2][1].j*hy;
    dx = fabs(wi[2][1].x-x);
    dy = fabs(wi[2][1].y-y);
    wi[2][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[2][2].i = rx + 1; wi[2][2].x = wi[2][2].i*hx;
    wi[2][2].j = ry + 1; wi[2][2].y = wi[2][2].j*hy;
    dx = fabs(wi[2][2].x-x);
    dy = fabs(wi[2][2].y-y);
    wi[2][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[1][2].i = rx + 1; wi[1][2].x = wi[1][2].i*hx;
    wi[1][2].j = ry + 0; wi[1][2].y = wi[1][2].j*hy;
    dx = fabs(wi[1][2].x-x);
    dy = fabs(wi[1][2].y-y);
    wi[1][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[0][0].i = rx - 1; wi[0][0].x = wi[0][0].i*hx;
    wi[0][0].j = ry - 1; wi[0][0].y = wi[0][0].j*hy;
    dx = fabs(wi[0][0].x-x);
    dy = fabs(wi[0][0].y-y);
    wi[0][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[1][0].i = rx - 1; wi[1][0].x = wi[1][0].i*hx;
    wi[1][0].j = ry + 0; wi[1][0].y = wi[1][0].j*hy;
    dx = fabs(wi[1][0].x-x);
    dy = fabs(wi[1][0].y-y);
    wi[1][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[2][0].i = rx - 1; wi[2][0].x = wi[2][0].i*hx;
    wi[2][0].j = ry + 1; wi[2][0].y = wi[2][0].j*hy;
    dx = fabs(wi[2][0].x-x);
    dy = fabs(wi[2][0].y-y);
    wi[2][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[3][0].i = rx - 1; wi[3][0].x = wi[3][0].i*hx;
    wi[3][0].j = ry + 2; wi[3][0].y = wi[3][0].j*hy;
    dx = fabs(wi[3][0].x-x);
    dy = fabs(wi[3][0].y-y);
    wi[3][0].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[3][1].i = rx + 0; wi[3][1].x = wi[3][1].i*hx;
    wi[3][1].j = ry + 2; wi[3][1].y = wi[3][1].j*hy;
    dx = fabs(wi[3][1].x-x);
    dy = fabs(wi[3][1].y-y);
    wi[3][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[3][2].i = rx + 1; wi[3][2].x = wi[3][2].i*hx;
    wi[3][2].j = ry + 2; wi[3][2].y = wi[3][2].j*hy;
    dx = fabs(wi[3][2].x-x);
    dy = fabs(wi[3][2].y-y);
    wi[3][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[3][3].i = rx + 2; wi[3][3].x = wi[3][3].i*hx;
    wi[3][3].j = ry + 2; wi[3][3].y = wi[3][3].j*hy;
    dx = fabs(wi[3][3].x-x);
    dy = fabs(wi[3][3].y-y);
    wi[3][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[2][3].i = rx + 2; wi[2][3].x = wi[2][3].i*hx;
    wi[2][3].j = ry + 1; wi[2][3].y = wi[2][3].j*hy;
    dx = fabs(wi[2][3].x-x);
    dy = fabs(wi[2][3].y-y);
    wi[2][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[1][3].i = rx + 2; wi[1][3].x = wi[1][3].i*hx;
    wi[1][3].j = ry + 0; wi[1][3].y = wi[1][3].j*hy;
    dx = fabs(wi[1][3].x-x);
    dy = fabs(wi[1][3].y-y);
    wi[1][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);

    wi[0][3].i = rx + 2; wi[0][3].x = wi[0][3].i*hx;
    wi[0][3].j = ry - 1; wi[0][3].y = wi[0][3].j*hy;
    dx = fabs(wi[0][3].x-x);
    dy = fabs(wi[0][3].y-y);
    wi[0][3].w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[0][2].i = rx + 1; wi[0][2].x = wi[0][2].i*hx;
    wi[0][2].j = ry - 1; wi[0][2].y = wi[0][2].j*hy;
    dx = fabs(wi[0][2].x-x);
    dy = fabs(wi[0][2].y-y);
    wi[0][2].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

    wi[0][1].i = rx + 0; wi[0][1].x = wi[0][1].i*hx;
    wi[0][1].j = ry - 1; wi[0][1].y = wi[0][1].j*hy;
    dx = fabs(wi[0][1].x-x);
    dy = fabs(wi[0][1].y-y);
    wi[0][1].w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
}

void ExtendedSpaceNode2D::clearWeights()
{
    for (unsigned int j=0; j<rows; j++) delete [] wi[j];
    delete [] wi;
}

void ExtendedSpaceNode2D::extendLayers(unsigned int layerNumber)
{
    this->layerNumber = layerNumber;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            wi[j][i].u = new double[layerNumber];
        }
    }
}

void ExtendedSpaceNode2D::clearLayers()
{
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            delete [] wi[j][i].u;
        }
    }
}

double ExtendedSpaceNode2D::value(unsigned int layer) const
{
    double P = 0.0;

    double Li = 0.0;
    double Lj = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

            P += Lj*Li*wi[j][i].u[layer];
        }
    }
    return P;
}

double ExtendedSpaceNode2D::value(double x, double y, unsigned int layer) const
{
    double P = 0.0;

    double Li = 0.0;
    double Lj = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

            P += Lj*Li*wi[j][i].u[layer];
        }
    }
    return P;
}

double ExtendedSpaceNode2D::valueDx(double x, double y, unsigned int layer) const
{
    double P = 0.0;

    double Li = 0.0;
    double Lj = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

            P += Lj*Li*wi[j][i].u[layer];
        }
    }
    return P;
}

double ExtendedSpaceNode2D::valueDy(double x, double y, unsigned int layer) const
{
    double P = 0.0;

    double Li = 0.0;
    double Lj = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

            P += Lj*Li*wi[j][i].u[layer];
        }
    }
    return P;
}

double ExtendedSpaceNode2D::valueDxN(unsigned int layer, double h) const
{
    return (value(x+h,y,layer)-value(x-h,y,layer))/(2.0*h);
}

double ExtendedSpaceNode2D::valueDyN(unsigned int layer, double h) const
{
    return (value(x,y+h,layer)-value(x,y-h,layer))/(2.0*h);
}

void P2Setting::toVector(DoubleVector &prms) const
{
    prms.clear();
    prms.resize(2*Lc*Lo + 2*Lo + 2*Lc);

    for (unsigned int i=0; i<Lc; i++)
    {
        prms[2*i + 0 + 2*Lc*Lo] = eta[i].x;
        prms[2*i + 1 + 2*Lc*Lo] = eta[i].y;

        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo + j] = k[i][j];
            prms[i*Lo + j + Lc*Lo] = z[i][j];

            prms[2*j + 0 + 2*Lc + 2*Lc*Lo] = xi[j].x;
            prms[2*j + 1 + 2*Lc + 2*Lc*Lo] = xi[j].y;
        }
    }
}

void P2Setting::fromVector(const DoubleVector &prms)
{
    k.clear();
    k.resize(Lc, Lo);

    z.clear();
    z.resize(Lc, Lo);

    xi.clear();
    xi.resize(Lo);

    eta.clear();
    eta.resize(Lc);

    for (unsigned int i=0; i<Lc; i++)
    {
        eta[i].x = prms[2*Lc*Lo + 2*i+0];
        eta[i].y = prms[2*Lc*Lo + 2*i+1];

        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = prms[i*Lo + j];
            z[i][j] = prms[i*Lo + j + Lc*Lo];

            xi[j].x = prms[2*Lc*Lo + 2*Lc + 2*j+0];
            xi[j].y = prms[2*Lc*Lo + 2*Lc + 2*j+1];
        }
    }
}

