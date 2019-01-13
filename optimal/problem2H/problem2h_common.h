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
#include "problem2h_global.h"

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
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;

#ifdef TIME_DISCRETE_H
    std::vector<double> tau;
#endif
};

struct InitialPulse
{
    double q;
    SpacePoint theta;
};

struct PROBLEM2HSHARED_EXPORT EquationParameterH
{
    double a;
    double lambda;

    unsigned int No;
    unsigned int Nc;

    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;

    unsigned int Ns;
    std::vector<SpacePoint> theta;
    std::vector<double> q;

//    std::vector<InitialPulse> pulseVector;

#ifdef TIME_DISCRETE_H
    unsigned int Nt;
#endif
};

struct PROBLEM2HSHARED_EXPORT EquationParameterH1
{
    struct InitialPulse
    {
        double q;
        SpacePoint theta;
    };

    struct OptimizationParam
    {
        DoubleMatrix k;
        DoubleMatrix z;
        std::vector<SpacePoint> xi;
        std::vector<SpacePoint> eta;
    };

    EquationParameterH1(unsigned int Nc, unsigned int No, unsigned int Ns, double a = 1.0, double alpha = 0.01);

    double a;
    double alpha;
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
private:
    OptimizationParam op;
    OptimizationParam rp;
    std::vector<InitialPulse> pulseVector;

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
