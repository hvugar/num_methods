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
};

struct PROBLEM2HSHARED_EXPORT SpacePointInfo : public SpaceNodePDE
{
    unsigned int id;
    double *u;
    double *ux;
    double *uy;
    unsigned int layerNumber;

    inline void clearWeights()
    {
        delete [] u;   u = NULL;
        delete [] ux;  ux = NULL;
        delete [] uy;  uy = NULL;
    }

    inline double value(unsigned int layer) const { return u[layer]; }
    inline double valueDx(unsigned int layer) const { return ux[layer]; }
    inline double valueDy(unsigned int layer) const { return uy[layer]; }
};

typedef std::vector<unsigned int>           uint_vector;
typedef std::set<unsigned int>              uint_set;
typedef std::vector<ExtendedSpacePointNode> espn_vector;
typedef std::vector<SpacePointInfo>         spif_vector;

#endif // PROBLEM2H_COMMON_H
