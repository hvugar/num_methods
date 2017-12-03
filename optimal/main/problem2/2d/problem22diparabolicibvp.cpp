#include "problem22diparabolicibvp.h"

double Problem22DIParabolicIBVP::boundary(const SpaceNodePDE&, const TimeNodePDE&, BoundaryType) const
{
    return NAN;
}

double Problem22DIParabolicIBVP::f(const SpaceNodePDE&, const TimeNodePDE&) const
{
    return 0.0;
}

void Problem22DIParabolicIBVP::setParametr(const Parameters &p)
{
    mParameters = p;
}

const Problem22DIParabolicIBVP::Parameters& Problem22DIParabolicIBVP::parametrs() const
{
    return mParameters;
}

void Problem22DIParabolicIBVP::extendDeltaControlPointsToGrid()
{
    for (unsigned int i=0; i<mParameters.Lc; i++)
    {
        extendDeltaControlPointToGrid1(mParameters.eta[i], i);
    }
}

void Problem22DIParabolicIBVP::extendDeltaControlPointToGrid1(ControlPoint &eta, unsigned int index)
{
    Dimension dimX = spaceDimension(Dimension::DimensionX);
    Dimension dimY = spaceDimension(Dimension::DimensionY);

    unsigned int NX = dimX.sizeN();
    unsigned int NY = dimY.sizeN();

    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int rx = (unsigned int)(round(eta.x*NX));
    unsigned int ry = (unsigned int)(round(eta.y*NY));

    double factor = 1.0 / (hx*hy);

    ExSpaceNodePDE espn;
    if ( rx*hx <= eta.x && ry*hy <= eta.y ) // left bottom
    {
        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);
    }
    else if ( rx*hx <= eta.x && ry*hy >= eta.y ) // left top
    {
        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);
    }
    else if ( rx*hx >= eta.x && ry*hy >= eta.y ) // right top
    {
        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);
        return;
    }
    else if ( rx*hx >= eta.x && ry*hy <= eta.y ) // right bottom
    {
        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        espn.w = ((hx-fabs(espn.x-eta.x))/hx)*((hy-fabs(espn.y-eta.y))/hy) * factor;
        eta.extNodes.push_back(espn);
    }
}

void Problem22DIParabolicIBVP::extendObservationPointToGrid()
{
    for (unsigned int j=0; j<mParameters.Lo; j++)
    {
        extendObservationPointToGrid1(mParameters.xi[j], j);
    }
}

void Problem22DIParabolicIBVP::extendObservationPointToGrid1(ObservationPoint &xi, unsigned int index)
{
    Dimension dimX = spaceDimension(Dimension::DimensionX);
    Dimension dimY = spaceDimension(Dimension::DimensionY);

    unsigned int NX = dimX.sizeN();
    unsigned int NY = dimY.sizeN();

    double hx = dimX.step();
    double hy = dimY.step();

    unsigned int rx = (unsigned int)(round(xi.x*NX));
    unsigned int ry = (unsigned int)(round(xi.y*NY));

    double hx30 = hx*hx*hx;
    double hx32 = 1.0 / (2.0*hx30);
    double hx36 = 1.0 / (6.0*hx30);

    double hy30 = hy*hy*hy;
    double hy32 = 1.0 / (2.0*hy30);
    double hy36 = 1.0 / (6.0*hy30);

    double dx = 0.0;
    double dy = 0.0;

    ExSpaceNodePDE espn;
    if ( rx*hx <= xi.x && ry*hy <= xi.y ) // left bottom
    {
        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36) * ((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx - 1; espn.x = espn.i*hx;
        espn.j = ry + 2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry + 2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry + 2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 2; espn.x = espn.i*hx;
        espn.j = ry + 2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 2; espn.x = espn.i*hx;
        espn.j = ry + 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 2; espn.x = espn.i*hx;
        espn.j = ry + 0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 2; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 1; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx + 0; espn.x = espn.i*hx;
        espn.j = ry - 1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);
    }
    else if ( rx*hx <= xi.x && ry*hy >= xi.y ) // left top
    {
        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);;
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+2; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+2; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+2; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+2; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);
    }
    else if ( rx*hx >= xi.x && ry*hy >= xi.y ) // right top
    {
        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);;
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry-2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);
    }
    else if ( rx*hx >= xi.x && ry*hy <= xi.y ) // right bottom
    {
        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-2; espn.x = espn.i*hx;
        espn.j = ry+2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry+2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);;
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry+2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+2; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry+0; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(hy+dy)*hy32);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(3.0*hx-dx)*hx36)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx+0; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);

        espn.index = index;
        espn.i = rx-1; espn.x = espn.i*hx;
        espn.j = ry-1; espn.y = espn.j*hy;
        dx = fabs(espn.x-xi.x);
        dy = fabs(espn.y-xi.y);
        espn.w = ((2.0*hx-dx)*(hx-dx)*(hx+dx)*hx32)*((2.0*hy-dy)*(hy-dy)*(3.0*hy-dy)*hy36);
        xi.extNodes.push_back(espn);
    }
}

void Problem22DIParabolicIBVP::Parameters::toVector(DoubleVector &prms) const
{
    prms.clear();
    prms.resize(2*Lc*Lo + 2*Lo + 2*Lc);

    // k
    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo + j] = k[i][j];
        }
    }

    // z
    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            prms[i*Lo + j + Lc*Lo] = z[i][j];
        }
    }

    // xi
    for (unsigned int j=0; j<Lo; j++)
    {
        prms[2*j + 0 + 2*Lc*Lo] = xi[j].x;
        prms[2*j + 1 + 2*Lc*Lo] = xi[j].y;
    }

    // eta
    for (unsigned int i=0; i<Lc; i++)
    {
        prms[2*i + 0 + 2*Lo + 2*Lc*Lo] = eta[i].x;
        prms[2*i + 1 + 2*Lo + 2*Lc*Lo] = eta[i].y;
    }
}

void Problem22DIParabolicIBVP::Parameters::fromVector(const DoubleVector &prms, unsigned int Lc, unsigned int Lo)
{
    this->Lc = Lc;
    this->Lo = Lo;

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
        for (unsigned int j=0; j<Lo; j++)
        {
            k[i][j] = prms[i*Lo + j];
        }
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        for (unsigned int j=0; j<Lo; j++)
        {
            z[i][j] = prms[Lc*Lo + i*Lo + j];
        }
    }

    for (unsigned int j=0; j<Lo; j++)
    {
        xi[j].x = prms[2*Lc*Lo + 2*j];
        xi[j].y = prms[2*Lc*Lo + 2*j+1];
    }

    for (unsigned int i=0; i<Lc; i++)
    {
        eta[i].x = prms[2*Lc*Lo + 2*Lo + 2*i];
        eta[i].y = prms[2*Lc*Lo + 2*Lo + 2*i+1];
    }
}
