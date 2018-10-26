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

struct PROBLEM2HSHARED_EXPORT OptimizeParameterH
{
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
};

struct PROBLEM2HSHARED_EXPORT EquationParameterH
{
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
    std::vector<SpacePoint> theta;
    DoubleVector q;
    double a;
    double lambda;
};

struct PROBLEM2HSHARED_EXPORT ExtendedSpacePointNodeH
{
    double x;
    double y;
    int nx;
    int ny;
    double w;
    bool isCenter;

    inline auto equals(const SpaceNodePDE &sn) const -> bool { return (sn.i == nx && sn.j == ny); }
};

struct PROBLEM2HSHARED_EXPORT ExtendedSpacePointH : public SpacePoint
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

    std::vector<ExtendedSpacePointNodeH> nodes;

    auto contains(int nx, int ny) const -> bool;
    auto contains(const SpaceNodePDE &sn) const -> bool { return ((minX <= sn.i && sn.i <= maxX) && (minY <= sn.j && sn.j <= maxY)); }
};

struct PROBLEM2HSHARED_EXPORT SpacePointInfoH: public SpacePoint
{
    SpacePointInfoH();
    SpacePointInfoH(unsigned int length);
    virtual ~SpacePointInfoH();
    void init(unsigned int length);
    void clear();

    unsigned int length;
    std::vector<double> vl;
    std::vector<double> dx;
    std::vector<double> dy;
};

typedef std::vector<unsigned int>            uint_vector;
typedef std::vector<SpacePointInfoH>         spif_vector;

#endif // PROBLEM2H_COMMON_H
