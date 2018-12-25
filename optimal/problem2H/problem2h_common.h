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
#include "problem2h_global.h"

//struct PROBLEM2HSHARED_EXPORT EquationParameterHE
//{
//    struct SpacePointExt
//    {
//        unsigned int rx;
//        unsigned int ry;
//        unsigned int minX;
//        unsigned int maxX;
//        unsigned int minY;
//        unsigned int maxY;

//        auto contains(int nx, int ny) const -> bool;
//        auto contains(const SpaceNodePDE &sn) const -> bool;

//        struct Node {
//            double x;
//            double y;
//            unsigned int nx;
//            unsigned int ny;
//            double w;
//            bool isCenter;
//        };

//        std::vector<Node> nodes;

//        unsigned int length;
//        std::vector<double> xi_vl;
//        std::vector<double> xi_dx;
//        std::vector<double> xi_dy;
//    };

//    EquationParameterHE();
//    EquationParameterHE(unsigned int Nc, unsigned int No, unsigned int Nd, unsigned int length, unsigned int gw = 3);

//    unsigned int Nc;
//    unsigned int No;
//    unsigned int Nd;
//#ifdef TIME_DISCRETE_H
//    unsigned int Nt;
//#endif
//    DoubleMatrix k;
//    DoubleMatrix z;
//    std::vector<SpacePoint> xi;
//    std::vector<SpacePoint> eta;
//#ifdef TIME_DISCRETE_H
//    DoubleVector tau;
//#endif
//    std::vector<SpacePoint> theta;
//    DoubleVector q;
//    double a;
//    double lambda;

//    std::vector<SpacePointExt> xi_ext;
//    std::vector<SpacePointExt> eta_ext;
//    std::vector<SpacePointExt> theta_ext;
//};

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

class DeltaGrid
{
public:
    DeltaGrid();
    virtual ~DeltaGrid();

    auto initGrid(unsigned int _N, double _hx, unsigned int _M, double _hy) -> void;
    auto setPoint(const SpacePoint &p, unsigned int sigmaX=1, unsigned int sigmaY=1) -> void;
    auto cleanGrid() -> void;

    auto isCenter(const SpaceNodePDE &sn) const -> bool;
    auto isCenter(unsigned int n, unsigned int m) const -> bool;

    auto isContains(const SpaceNodePDE &sn) const -> bool;
    auto isContains(unsigned int n, unsigned int m) const -> bool;

    auto weight(const SpaceNodePDE &sn) const -> double;
    auto weight(unsigned int n, unsigned int m) const -> double;

    unsigned int minX() const { return _minX; }
    unsigned int maxX() const { return _maxX; }
    unsigned int minY() const { return _minY; }
    unsigned int maxY() const { return _maxY; }

    double **data() { return nwmx; }
    double **data() const { return nwmx; }

private:

    unsigned int _N;
    unsigned int _M;
    double _hx;
    double _hy;
    SpacePoint _p;

    double **nwmx;
    unsigned int _rx;
    unsigned int _ry;
    unsigned int _minX;
    unsigned int _maxX;
    unsigned int _minY;
    unsigned int _maxY;
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

struct PROBLEM2HSHARED_EXPORT EquationParameterH
{
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
    std::vector<SpacePoint> theta;
    DoubleVector q;
    double a;
    double lambda;
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
