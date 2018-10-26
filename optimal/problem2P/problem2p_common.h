#ifndef PROBLEM2P_COMMON_H
#define PROBLEM2P_COMMON_H

#include "problem2p_global.h"
#include <vector2d.h>
#include <projection.h>
#include <printer.h>
#include <grid/ibvp.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <vectornormalizer.h>

#include "../imaging/imaging.h"

#include <algorithm>
#include <vector>
#include <utility>

struct OptimizeParameterP
{
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
};

struct EquationParameterP
{
    unsigned int No;
    unsigned int Nc;
    double a;
    double alpha;
    double lambda;
    double theta;
    double phi;
};

struct ExtendedSpacePointNodeP
{
    double x;
    double y;
    int nx;
    int ny;
    double w;
    bool isCenter;

    inline auto equals(const SpaceNodePDE &sn) const -> bool { return (sn.i == nx && sn.j == ny); }
};

struct ExtendedSpacePointP : public SpacePoint
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

    std::vector<ExtendedSpacePointNodeP> nodes;

    auto contains(int nx, int ny) const -> bool;
    auto contains(const SpaceNodePDE &sn) const -> bool
    {
        return ((minX <= sn.i && sn.i <= maxX) &&
                (minY <= sn.j && sn.j <= maxY));
    }
};

struct SpacePointInfoP : public SpacePoint
{
    SpacePointInfoP();
    SpacePointInfoP(unsigned int length);
    virtual ~SpacePointInfoP();
    void init(unsigned int length);
    void clear();

    unsigned int length;
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
};

typedef std::vector<unsigned int>     uint_vector;
typedef std::vector<SpacePointInfoP>  spif_vector;

#endif // PROBLEM2P_COMMON_H
