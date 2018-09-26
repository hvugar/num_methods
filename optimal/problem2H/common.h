#ifndef PROBLEM2H_COMMON_H
#define PROBLEM2H_COMMON_H

#include <vector>
#include <set>
#include <algorithm>
#include <time.h>
#include <function.h>
#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <grid/ibvp.h>
#include "problem2h_global.h"
#include "../imaging/imaging.h"

//double hx = dimX.step();
//double hy = dimY.step();

//unsigned int Nx = dimX.sizeN();
//unsigned int Ny = dimY.sizeN();

//unsigned int rx = (unsigned int) ( round(pt.x*Nx) );
//unsigned int ry = (unsigned int) ( round(pt.y*Ny) );

//double sigmaX = hx;
//double sigmaY = hy;

//double sumX = 0.0;
//for (unsigned int n=rx-k; n<=rx+k; n++) sumX += exp(-((n*hx-pt.x)*(n*hx-pt.x))/(2.0*sigmaX*sigmaX));
//sumX *= hx;

//double sumY = 0.0;
//for (unsigned int m=ry-k; m<=ry+k; m++) sumY += exp(-((m*hy-pt.y)*(m*hy-pt.y))/(2.0*sigmaY*sigmaY));
//sumY *= hy;

//double sigma = (sumX*sumY) / (2.0*M_PI);
//double factor = 1.0/((2.0*M_PI)*sigma);

//for (unsigned int m=ry-k; m<=ry+k; m++)
//{
//    for (unsigned int n=rx-k; n<=rx+k; n++)
//    {
//        ExtendedSpacePointNode node;
//        node.i = n; node.x = n*hx;
//        node.j = m; node.y = m*hy;
//        node.pt = pt; node.id = id;
//        node.w = factor*exp(-0.5*(((node.x-pt.x)*(node.x-pt.x))/(sigmaX*sigmaX)+((node.y-pt.y)*(node.y-pt.y))/(sigmaY*sigmaY)));
//        node.isCenter = ( m==ry && n==rx );
//        nodes.push_back(node);
//    }
//}


struct PulseInfluence
{
    unsigned int id = 0;
    double x;
    double y;
    unsigned int rx;
    unsigned int ry;
    double q;
    unsigned int extCount;
    double **w;

    bool contains(unsigned int n, unsigned int m) const { return (rx-extCount <= n && n <= rx+extCount) && (ry-extCount <= m && m <= ry+extCount);  }

    double q_w(unsigned int n, unsigned m) const { return contains(n,m) ? w[m][n] : 0.0; }

    void distribute(double x, double y, unsigned int Nx, unsigned int Ny/*, double hx, double hy*/)
    {
        rx = (unsigned int) round( x*Nx );
        ry = (unsigned int) round( y*Ny );

        w = new double*[2*extCount + 1];
        for (unsigned int m = 0; m <= 2*extCount; m++)
        {
            w[m] = new double[2*extCount + 1];
        }
    }
};

struct PROBLEM2HSHARED_EXPORT OptimizeParameter
{
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
};

struct PROBLEM2HSHARED_EXPORT EquationParameter
{
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
    std::vector<SpacePoint> theta;
    DoubleVector q;
    double a;
    double lambda;
};

struct PROBLEM2HSHARED_EXPORT ExtendedSpacePointNode
{
    SpacePoint pt;
    unsigned int id;
    unsigned int i;
    unsigned int j;
    double x;
    double y;
    double w;
    bool isCenter;
};

struct PROBLEM2HSHARED_EXPORT SpacePointInfo : public SpaceNodePDE
{
    SpacePointInfo();
    SpacePointInfo(unsigned int id, unsigned int layerNumber);
    virtual ~SpacePointInfo();

    void createSpacePointInfos(unsigned int layerNumber);
    void clearWeights();

    inline double vl(unsigned int layer) const { return _vl[layer]; }
    inline double dx(unsigned int layer) const { return _dx[layer]; }
    inline double dy(unsigned int layer) const { return _dy[layer]; }

    unsigned int id;
    unsigned int layerNumber;

    double *_vl;
    double *_dx;
    double *_dy;
};

typedef std::vector<unsigned int>           uint_vector;
typedef std::set<unsigned int>              uint_set;
typedef std::vector<ExtendedSpacePointNode> espn_vector;
typedef std::vector<SpacePointInfo>         spif_vector;

#endif // PROBLEM2H_COMMON_H
