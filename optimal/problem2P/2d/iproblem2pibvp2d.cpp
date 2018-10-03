#include "iproblem2pibvp2d.h"

void extendPointToGridNodes(const SpaceNodePDE &point, int id, vector<ExtendedGridNode> nodes, const Dimension &dimensionX, const Dimension &dimensionY)
{
    unsigned int Nx = dimensionX.sizeN();
    unsigned int Ny = dimensionY.sizeN();

    double hx = dimensionX.step();
    double hy = dimensionY.step();

    unsigned int rx = (unsigned int)(floor(point.x*Nx));
    unsigned int ry = (unsigned int)(floor(point.y*Ny));

    double hx3 = hx*hx*hx;
    double hx32 = (1.0/(2.0*hx3));
    double hx36 = (1.0/(6.0*hx3));

    double hy3 = hy*hy*hy;
    double hy32 = 1.0/(2.0*hy3);
    double hy36 = 1.0/(6.0*hy3);

    ExtendedGridNode node;
    double dx = 0.0;
    double dy = 0.0;

    if (rx == 0 && ry == 0)
    {}
    else if (rx == 0 && ry == Ny)
    {}
    else if (rx == Nx && ry == Ny)
    {}
    else if (rx == Nx && ry == 0)
    {}
    else if (rx==0 && ry == 1)
    {}
    else if (rx==1 && ry == 0)
    {}
    else if (rx == 0 && ry == Ny-1)
    {}
    else if (rx == 1 && ry == Ny)
    {}
    else
    {
        node.i = rx-1; node.x = node.i*hx;
        node.j = ry-1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);

        node.i = rx-1; node.x = node.i*hx;
        node.j = ry+0; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx-1; node.x = node.i*hx;
        node.j = ry+1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx-1; node.x = node.i*hx;
        node.j = ry+2; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+0; node.x = node.i*hx;
        node.j = ry-1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+0; node.x = node.i*hx;
        node.j = ry+0; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+0; node.x = node.i*hx;
        node.j = ry+1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+0; node.x = node.i*hx;
        node.j = ry+2; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+1; node.x = node.i*hx;
        node.j = ry-1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+1; node.x = node.i*hx;
        node.j = ry+0; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+1; node.x = node.i*hx;
        node.j = ry+1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+1; node.x = node.i*hx;
        node.j = ry+2; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+2; node.x = node.i*hx;
        node.j = ry-1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);

        node.i = rx+2; node.x = node.i*hx;
        node.j = ry+0; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+2; node.x = node.i*hx;
        node.j = ry+1; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        nodes.push_back(node);

        node.i = rx+2; node.x = node.i*hx;
        node.j = ry+2; node.y = node.j*hy;
        node.point = point; node.id = id;
        dx = fabs(node.x-point.x); dy = fabs(node.y-point.y); node.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        nodes.push_back(node);
    }
}

//-------------------------------------------------------------------------------------------------------//

double IProblem22DPIBVP::boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryType) const
{
    return NAN;
}

double IProblem22DPIBVP::f(const SpaceNodePDE&, const TimeNodePDE&) const
{
    return 0.0;
}

void IProblem22DPIBVP::setEquationParameters(double a, double lambda0, double lambda)
{
    this->a = a;
    this->lambda0 = lambda0;
    this->lambda = lambda;
}

void IProblem22DPIBVP::setIntTemperature(double fi)
{
    this->fi = fi;
}

void IProblem22DPIBVP::setEnvTemperature(double theta)
{
    this->theta = theta;
}

void IProblem22DPIBVP::setParameter(const Parameter &parameter)
{
    mParameter = parameter;
}

const Parameter& IProblem22DPIBVP::parameter() const
{
    return mParameter;
}

void IProblem22DPIBVP::distributeDelta(const SpaceNodePDE &pt, std::vector<ExtendedDeltaPoint> &nodes, unsigned int id) const
{
    Dimension xd = spaceDimension(Dimension::DimensionX);
    Dimension yd = spaceDimension(Dimension::DimensionY);
    double hx = xd.step();
    double hy = yd.step();
    unsigned int Nx = xd.sizeN();
    unsigned int Ny = yd.sizeN();

    double sigmaX = hx;
    double sigmaY = hy;

    unsigned int rx = (unsigned int)(round(pt.x*Nx));
    unsigned int ry = (unsigned int)(round(pt.y*Ny));

    unsigned int k=3;

    double sumX = 0.0;
    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
    }
    sumX *= hx;

    double sumY = 0.0;
    for (unsigned int m=ry-k; m<=ry+k; m++)
    {
        sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
    }
    sumY *= hy;

    double sigma = (sumX*sumY) / (2.0*M_PI);
    double factor = 1.0/((2.0*M_PI)*sigma);
    //double factor = (1.0/(2.0*M_PI*sigmaX*sigmaY));

    for (unsigned int n=rx-k; n<=rx+k; n++)
    {
        for (unsigned int m=ry-k; m<=ry+k; m++)
        {
            ExtendedDeltaPoint node;
            node.i = n; node.x = n*hx; node.j = m; node.y = m*hy; node.pt = pt; node.id = id;
            node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
            nodes.push_back(node);
        }
    }
}

//-------------------------------------------------------------------------------------------------------//

ExtendedSpaceNode2D::ExtendedSpaceNode2D()
{
    x = y = 0.0;
    i = j = 0;
    wi = NULL;
    layerNumber = 0;
}

ExtendedSpaceNode2D::~ExtendedSpaceNode2D()
{}

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
    for (unsigned int rw=0; rw<rows; rw++)
        wi[rw] = new WISpaceNodePDE[cols];

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
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double P = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            P += Ly[j]*Lx[i]*u;
        }
    }
    return P;
}

double ExtendedSpaceNode2D::valueDx(unsigned int layer) const
{
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)+(x-x2)*(x-x3)+(x-x1)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)+(x-x2)*(x-x3)+(x-x0)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)+(x-x1)*(x-x3)+(x-x0)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)+(x-x1)*(x-x2)+(x-x0)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double Px = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            Px += Ly[j]*Lx[i]*u;
        }
    }
    return Px;
}

double ExtendedSpaceNode2D::valueDy(unsigned int layer) const
{
    double Lx[] = {0.0, 0.0, 0.0, 0.0};
    double Ly[] = {0.0, 0.0, 0.0, 0.0};

    double y0 = wi[0][0].y;
    double y1 = wi[1][0].y;
    double y2 = wi[2][0].y;
    double y3 = wi[3][0].y;
    Ly[0] = ((y-y1)*(y-y2)+(y-y2)*(y-y3)+(y-y1)*(y-y3))/((y0-y1)*(y0-y2)*(y0-y3));
    Ly[1] = ((y-y0)*(y-y2)+(y-y2)*(y-y3)+(y-y0)*(y-y3))/((y1-y0)*(y1-y2)*(y1-y3));
    Ly[2] = ((y-y0)*(y-y1)+(y-y1)*(y-y3)+(y-y0)*(y-y3))/((y2-y0)*(y2-y1)*(y2-y3));
    Ly[3] = ((y-y0)*(y-y1)+(y-y1)*(y-y2)+(y-y0)*(y-y2))/((y3-y0)*(y3-y1)*(y3-y2));

    double x0 = wi[0][0].x;
    double x1 = wi[0][1].x;
    double x2 = wi[0][2].x;
    double x3 = wi[0][3].x;
    Lx[0] = ((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3));
    Lx[1] = ((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3));
    Lx[2] = ((x-x0)*(x-x1)*(x-x3))/((x2-x0)*(x2-x1)*(x2-x3));
    Lx[3] = ((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2));

    double Py = 0.0;
    for (unsigned int j=0; j<rows; j++)
    {
        for (unsigned int i=0; i<cols; i++)
        {
            double u = wi[j][i].u[layer];
            Py += Ly[j]*Lx[i]*u;
        }
    }
    return Py;
}

//double ExtendedSpaceNode2D::value(double x, double y, unsigned int layer) const
//{
//    double P = 0.0;

//    double Li = 0.0;
//    double Lj = 0.0;
//    for (unsigned int j=0; j<rows; j++)
//    {
//        for (unsigned int i=0; i<cols; i++)
//        {
//            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
//            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
//            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
//            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


//            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
//            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
//            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
//            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

//            P += Lj*Li*wi[j][i].u[layer];
//        }
//    }
//    return P;
//}

//double ExtendedSpaceNode2D::valueDx(double x, double y, unsigned int layer) const
//{
//    double P = 0.0;

//    double Li = 0.0;
//    double Lj = 0.0;
//    for (unsigned int j=0; j<rows; j++)
//    {
//        for (unsigned int i=0; i<cols; i++)
//        {
//            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
//            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
//            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
//            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


//            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
//            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
//            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
//            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

//            P += Lj*Li*wi[j][i].u[layer];
//        }
//    }
//    return P;
//}

//double ExtendedSpaceNode2D::valueDy(double x, double y, unsigned int layer) const
//{
//    double P = 0.0;

//    double Li = 0.0;
//    double Lj = 0.0;
//    for (unsigned int j=0; j<rows; j++)
//    {
//        for (unsigned int i=0; i<cols; i++)
//        {
//            if (j==0) Lj = (((y-wi[1][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[0][i].y-wi[1][i].y)*(wi[0][i].y-wi[2][i].y)*(wi[0][i].y-wi[3][i].y)));
//            if (j==1) Lj = (((y-wi[0][i].y)*(y-wi[2][i].y)*(y-wi[3][i].y))/((wi[1][i].y-wi[0][i].y)*(wi[1][i].y-wi[2][i].y)*(wi[1][i].y-wi[3][i].y)));
//            if (j==2) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[3][i].y))/((wi[2][i].y-wi[0][i].y)*(wi[2][i].y-wi[1][i].y)*(wi[2][i].y-wi[3][i].y)));
//            if (j==3) Lj = (((y-wi[0][i].y)*(y-wi[1][i].y)*(y-wi[2][i].y))/((wi[3][i].y-wi[0][i].y)*(wi[3][i].y-wi[1][i].y)*(wi[3][i].y-wi[2][i].y)));


//            if (i==0) Li = (((x-wi[j][1].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][0].x-wi[j][1].x)*(wi[j][0].x-wi[j][2].x)*(wi[j][0].x-wi[j][3].x)));
//            if (i==1) Li = (((x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][2].x)*(x-wi[j][3].x))/((wi[j][1].x-wi[j][0].x)*(wi[j][1].x-wi[j][2].x)*(wi[j][1].x-wi[j][3].x)));
//            if (i==2) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][3].x)+(x-wi[j][1].x)*(x-wi[j][3].x))/((wi[j][2].x-wi[j][0].x)*(wi[j][2].x-wi[j][1].x)*(wi[j][2].x-wi[j][3].x)));
//            if (i==3) Li = (((x-wi[j][0].x)*(x-wi[j][1].x)+(x-wi[j][0].x)*(x-wi[j][2].x)+(x-wi[j][1].x)*(x-wi[j][2].x))/((wi[j][3].x-wi[j][0].x)*(wi[j][3].x-wi[j][1].x)*(wi[j][3].x-wi[j][2].x)));

//            P += Lj*Li*wi[j][i].u[layer];
//        }
//    }
//    return P;
//}

//double ExtendedSpaceNode2D::valueDxN(unsigned int layer, double h) const
//{
//    return (value(x+h,y,layer)-value(x-h,y,layer))/(2.0*h);
//}

//double ExtendedSpaceNode2D::valueDyN(unsigned int layer, double h) const
//{
//    return (value(x,y+h,layer)-value(x,y-h,layer))/(2.0*h);
//}

Parameter::Parameter(unsigned int Lc, unsigned int Lo) : Lc(Lc), Lo(Lo)
{
    k.resize(Lc, Lo);
    z.resize(Lc, Lo);
    eta.resize(Lc);
    xi.resize(Lo);
}

Parameter::Parameter(const DoubleVector &prmtrs, unsigned int Lc, unsigned int Lo) : Lc(Lc), Lo(Lo)
{
    k.resize(Lc, Lo);
    z.resize(Lc, Lo);
    eta.resize(Lc);
    xi.resize(Lo);

    for (unsigned int i=0; i<Lc; i++)
    {
        eta[i].x = prmtrs[2*Lc*Lo + 2*i+0];
        eta[i].y = prmtrs[2*Lc*Lo + 2*i+1];

        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = prmtrs[i*Lo + j];
            z[i][j] = prmtrs[i*Lo + j + Lc*Lo];

            xi[j].x = prmtrs[2*Lc*Lo + 2*Lc + 2*j+0];
            xi[j].y = prmtrs[2*Lc*Lo + 2*Lc + 2*j+1];
        }
    }
}

//void Parameter::toVector(DoubleVector &prms) const
//{
//    prms.clear();
//    prms.resize(2*Lc*Lo + 2*Lc + 2*Lo);

//    for (unsigned int i=0; i<Lc; i++)
//    {
//        prms[2*i + 0 + 2*Lc*Lo] = eta[i].x;
//        prms[2*i + 1 + 2*Lc*Lo] = eta[i].y;

//        for (unsigned int j=0; j<Lo; j++)
//        {
//            prms[i*Lo + j] = k[i][j];
//            prms[i*Lo + j + Lc*Lo] = z[i][j];

//            prms[2*j + 0 + 2*Lc + 2*Lc*Lo] = xi[j].x;
//            prms[2*j + 1 + 2*Lc + 2*Lc*Lo] = xi[j].y;
//        }
//    }
//}

//void Parameter::fromVector(const DoubleVector &prms)
//{
//    k.clear();
//    k.resize(Lc, Lo);

//    z.clear();
//    z.resize(Lc, Lo);

//    eta.clear();
//    eta.resize(Lc);

//    xi.clear();
//    xi.resize(Lo);

//    for (unsigned int i=0; i<Lc; i++)
//    {
//        eta[i].x = prms[2*Lc*Lo + 2*i+0];
//        eta[i].y = prms[2*Lc*Lo + 2*i+1];

//        for (unsigned int j=0; j<Lo; j++)
//        {
//            k[i][j] = prms[i*Lo + j];
//            z[i][j] = prms[i*Lo + j + Lc*Lo];

//            xi[j].x = prms[2*Lc*Lo + 2*Lc + 2*j+0];
//            xi[j].y = prms[2*Lc*Lo + 2*Lc + 2*j+1];
//        }
//    }
//}
