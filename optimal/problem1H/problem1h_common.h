#ifndef PROBLEM1H_COMMON_H
#define PROBLEM1H_COMMON_H

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
#include <benchmark.h>
#include <grid/ibvp.h>
#include <grid/hibvp.h>
#include <deltagrid.h>
#include "problem1h_global.h"

struct PROBLEM1HSHARED_EXPORT InitialPulse1D
{
    InitialPulse1D();
    InitialPulse1D(const SpacePoint &sp, double q);

    SpacePoint theta;
    double q;
};

struct PROBLEM1HSHARED_EXPORT OptimizeParameter1H
{
#if defined(DISCRETE_DELTA_TIME)
    DoubleMatrix *k;
    DoubleMatrix *z;
#else
    DoubleMatrix k;
    DoubleMatrix z;
#endif
    std::vector<SpacePoint> ksi;
    std::vector<SpacePoint> eta;
};

struct PROBLEM1HSHARED_EXPORT EquationParameter1H
{
    double a;
    double alpha;

    unsigned int No;
    unsigned int Nc;

#if defined(DISCRETE_DELTA_TIME)
    unsigned int Nt;
    std::vector<double> timeMoments;
#endif

    unsigned int Ns;
    std::vector<InitialPulse1D> pulses;

    DoubleVector Q1;
    DoubleVector Q2;
    DoubleVector X1;
    DoubleVector X2;

#ifdef TIME_DISCRETE_H
    unsigned int Nt;
#endif
};

struct PROBLEM1HSHARED_EXPORT ExtendedSpacePointNode1H
{
    double x;
    int nx;
    double w;
    bool isCenter;

    inline auto equals(const SpaceNodePDE &sn) const -> bool { return (sn.i == nx); }
};

struct PROBLEM1HSHARED_EXPORT ExtendedSpacePoint1H : public SpacePoint
{
    unsigned int id;
    int rx;

    int k;
    int minX;
    int maxX;

    std::vector<ExtendedSpacePointNode1H> nodes;

    auto contains(int nx, int ny) const -> bool;
    auto contains(const SpaceNodePDE &sn) const -> bool { return ((minX <= sn.i && sn.i <= maxX)); }
};

struct PROBLEM1HSHARED_EXPORT SpacePointInfo1H: public SpacePoint
{
    SpacePointInfo1H();
    SpacePointInfo1H(unsigned int length);
    virtual ~SpacePointInfo1H();
    void init(unsigned int length);
    void clear();

    unsigned int length;
    std::vector<double> vl;
    std::vector<double> dx;
};

typedef std::vector<unsigned int>     uint_vector1H;
typedef std::vector<SpacePointInfo1H> spif_vector1H;

#endif // PROBLEM2H_COMMON_H
