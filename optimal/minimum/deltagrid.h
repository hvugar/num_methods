#ifndef DELTA_GRID_H
#define DELTA_GRID_H

#include "grid/grid.h"
#include <cmath>
#include <exception>

class MINIMUMSHARED_EXPORT DeltaGridException : public std::exception
{
public:
    virtual const char* what() const _NOEXCEPT;
};

class MINIMUMSHARED_EXPORT DeltaGrid2D
{
public:
    DeltaGrid2D();
    virtual ~DeltaGrid2D();

    auto initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void;
    auto cleanGrid() -> void;

    auto distributeGauss(const SpacePoint& sp, unsigned sigmaXNum = 1, unsigned int sigmaYNum = 1) -> void;
    auto distributeSigle(const SpacePoint& sp) -> void;
    auto distributeRect4(const SpacePoint& sp) -> void;

    auto isCenter(const SpaceNodePDE &sn) const -> bool;
    auto isCenter(unsigned int n, unsigned int m) const -> bool;

    auto isContains(const SpaceNodePDE &sn) const -> bool;
    auto isContains(unsigned int n, unsigned int m) const -> bool;

    auto weight(const SpaceNodePDE &sn) const -> double;
    auto weight(unsigned int n, unsigned int m) const -> double;

    auto rx() const -> unsigned int;
    auto ry() const -> unsigned int;

    auto p() const -> const SpacePoint&;
    auto p() -> SpacePoint&;

    auto minX() const -> unsigned int;
    auto maxX() const -> unsigned int;
    auto minY() const -> unsigned int;
    auto maxY() const -> unsigned int;

    double **nodes() { return m_nodes; }
    double **nodes() const { return m_nodes; }

private:

    unsigned int mN;
    unsigned int mM;
    double mhx;
    double mhy;
    SpacePoint mp;

    double **m_nodes;
    unsigned int _rx;
    unsigned int _ry;
    unsigned int _minX;
    unsigned int _maxX;
    unsigned int _minY;
    unsigned int _maxY;
};

#endif // DELTA_GRID_H
