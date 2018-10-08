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

struct PROBLEM2HSHARED_EXPORT ExtendedSpacePointNode1
{
    double x;
    double y;
    int nx;
    int ny;
    double w;
    bool isCenter;

    inline auto equals(const SpaceNodePDE &sn) const -> bool { return (sn.i == nx && sn.j == ny); }
};

struct PROBLEM2HSHARED_EXPORT ExtendedSpacePoint : public SpacePoint
{
    unsigned int id;
    int rx;
    int ry;
    int rz;

    int k;
    int minX;
    int maxX;
    int minY;
    int maxY;

    std::vector<ExtendedSpacePointNode1> nodes;

    auto contains(int nx, int ny) const -> bool;
    auto contains(const SpaceNodePDE &sn) const -> bool { return ((minX <= sn.i && sn.i <= maxX) && (minY <= sn.j && sn.j <= maxY)); }
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

struct PROBLEM2HSHARED_EXPORT SpacePointInfo : public SpacePoint
{
    SpacePointInfo();
    SpacePointInfo(unsigned int length);
    virtual ~SpacePointInfo();
    void init(unsigned int length);
    void clear();

    unsigned int length;
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
};

//struct SpacePointExt : public SpacePoint
//{
//    typedef std::pair<unsigned int, unsigned int> GridNodePair;

//    struct GridNodeWeight
//    {
//        unsigned int id;
//        unsigned int i;
//        unsigned int j;
//        double x;
//        double y;
//        double w;
//        bool isCenter;
//    };
//    typedef std::pair<GridNodePair, GridNodeWeight> GridNodePairWeight;
//    typedef std::map <GridNodePair, GridNodeWeight> GridNodeMap;

//    unsigned int id;
//    unsigned int rx;
//    unsigned int ry;

//    const GridNodeMap &distPoints() const { return m_distPoints; }
//    GridNodeMap &distPoints() { return m_distPoints; }

//private:
//    GridNodeMap m_distPoints;
//};

//typedef std::vector<SpacePointExt> vector_SpacePointExt;
typedef std::vector<unsigned int>           uint_vector;
typedef std::set<unsigned int>              uint_set;
typedef std::vector<ExtendedSpacePointNode> espn_vector;
typedef std::vector<SpacePointInfo>         spif_vector;

#endif // PROBLEM2H_COMMON_H
