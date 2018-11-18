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
#include <benchmark.h>

#include "../imaging/imaging.h"

#include <algorithm>
#include <vector>
#include <utility>

#define TIME_DISCRETE1
#define OH1

struct OptimizeParameterP
{
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
#ifdef TIME_DISCRETE
    std::vector<double> tau;
#endif
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
#ifdef TIME_DISCRETE
    unsigned int Nt;
#endif
};

struct EquationParameterPr
{
    EquationParameterPr(double a, double alpha, double lambda, double theta, double phi, unsigned int Nc, unsigned int No)
        : a(a), alpha(alpha), lambda(lambda), theta(theta), phi(phi), Nc(Nc), No(No) {}

    double a;
    double alpha;
    double lambda;
    double theta;
    double phi;

    unsigned int Nc = 0;
    unsigned int No = 0;

    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;

#ifdef TIME_DISCRETE
    unsigned int Nt;
    std::vector<double> tau;
#endif

    auto toVector(DoubleVector &x) const -> void;
    auto fromVector(const DoubleVector &x, unsigned int Nc, unsigned int No) -> void;
};

struct ExtendedSpacePointNodeP
{
    double x;
    double y;
    int nx;
    int ny;
    double w;
    bool isCenter;

    inline auto equals(const SpaceNodePDE &sn) const -> bool { return (sn.i == static_cast<unsigned int>(nx) && sn.j == static_cast<unsigned int>(ny)); }
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
        return ((minX <= static_cast<int>(sn.i) && static_cast<int>(sn.i) <= maxX) &&
                (minY <= static_cast<int>(sn.j) && static_cast<int>(sn.j) <= maxY));
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
    std::vector<double> dxx;
    std::vector<double> dyy;
};

typedef std::vector<unsigned int>     uint_vector;
typedef std::vector<SpacePointInfoP>  spif_vector;

class grid_exception : public std::exception
{
public:
    grid_exception(const std::string &message) : std::exception() { this->message = message; }
    virtual char const* what() const _GLIBCXX_USE_NOEXCEPT { return message.data(); }
private:
    std::string message;
};

#endif // PROBLEM2P_COMMON_H
