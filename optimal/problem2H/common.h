#ifndef PROBLEM2H_COMMON_H
#define PROBLEM2H_COMMON_H

#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <utility>
#include <time.h>
#include <function.h>
#include <gradient.h>
#include <projection.h>
#include <printer.h>
#include <gradient_cjt.h>
#include <gradient_sd.h>
#include <gradient_cs.h>
#include <vector2d.h>
#include <vectornormalizer.h>
#include <grid/ibvp.h>
#include "problem2h_global.h"
#include "../imaging/imaging.h"

struct SpacePointExt : public SpacePoint
{
    typedef std::pair<unsigned int, unsigned int> GridNodePair;

    struct GridNodeWeight
    {
        unsigned int id;
        unsigned int i;
        unsigned int j;
        double x;
        double y;
        double w;
        bool isCenter;
    };
    typedef std::pair<GridNodePair, GridNodeWeight> GridNodePairWeight;
    typedef std::map <GridNodePair, GridNodeWeight> GridNodeMap;

    unsigned int id;
    unsigned int rx;
    unsigned int ry;
    GridNodeMap distPoints;
};

typedef std::vector<SpacePointExt> vector_SpacePointExt;

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
