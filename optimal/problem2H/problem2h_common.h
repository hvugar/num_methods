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
#include <benchmark.h>
#include <grid/ibvp.h>
#include <grid/hibvp.h>
#include <deltagrid.h>
#include <utils/random.h>
#include "problem2h_global.h"

struct PROBLEM2HSHARED_EXPORT InitialPulse2D
{
    InitialPulse2D();
    InitialPulse2D(const SpacePoint &sp, double q);
    InitialPulse2D(const InitialPulse2D& pulse);

    SpacePoint theta;
    double q;
};

struct PROBLEM2HSHARED_EXPORT EquaParameter2H
{
    double a;
    double alpha;

    unsigned int Np;
    std::vector<InitialPulse2D> pulses;
    unsigned int Nt;
    std::vector<TimeNodePDE> tm;

    unsigned int Nc;
    unsigned int No;
    struct {
        std::vector<DoubleMatrix> k;
        std::vector<DoubleMatrix> z;
        std::vector<SpacePoint> ksi;
        std::vector<SpacePoint> eta;
    } opt, reg;

    auto initPulseParemeters(unsigned int Np) -> void;
    auto initParemeters(unsigned int Nt, unsigned int Nc, unsigned int No) -> void;
    auto clearParameters() -> void;

    auto OptimalParameterFromVector(const DoubleVector &x) -> void;
    auto OptimalParameterToVector(DoubleVector &x) const -> void;

    auto RegularParameterToVector(DoubleVector &x) const -> void;
    auto RegularParameterFromVector(const DoubleVector &x) -> void;

    auto printOptimalParemeters() const -> void;
    auto printRegularParemeters() const -> void;
};

struct PROBLEM2HSHARED_EXPORT FuncParameter2H
{
    DoubleVector vmin;
    DoubleVector vmax;
    double r = 0.0;
    double regEpsilon = 0.0;
    DoubleVector Q1;
    DoubleVector Q2;
    std::vector<SpacePoint> Theta1;
    std::vector<SpacePoint> Theta2;
};

struct GridH
{
    std::vector<unsigned int> rows0;
    std::vector<unsigned int> rows1;
    std::vector<unsigned int> rows2;

    std::vector<unsigned int> cols0;
    std::vector<unsigned int> cols1;
    std::vector<unsigned int> cols2;
};

struct GridSpace2D
{
    unsigned int N;
    unsigned int M;
    double hx;
    double hy;
};

struct PROBLEM2HSHARED_EXPORT OptimizeParameterH
{
#if defined (DISCRETE_DELTA_TIME)
    DoubleMatrix *k;
    DoubleMatrix *z;
#else
    DoubleMatrix k;
    DoubleMatrix z;
#endif
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;

#ifdef TIME_DISCRETE_H
    std::vector<double> tau;
#endif
};

struct PROBLEM2HSHARED_EXPORT EquationParameterH
{
    double a;
    double alpha;

    unsigned int No;
    unsigned int Nc;

#if defined (DISCRETE_DELTA_TIME)
    unsigned int Nt;
    std::vector<TimeNodePDE> timeMoments;
#endif

    unsigned int Ns;
    std::vector<InitialPulse2D> pulses;

    DoubleVector Q1;
    DoubleVector Q2;

#ifdef TIME_DISCRETE_H
    unsigned int Nt;
#endif
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

typedef std::vector<unsigned int>    uint_vectorH;
typedef std::vector<SpacePointInfoH> spif_vectorH;

#endif // PROBLEM2H_COMMON_H
