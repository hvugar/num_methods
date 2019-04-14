#ifndef DELTA_GRID_H
#define DELTA_GRID_H

#include "grid/grid.h"
#include <cmath>
#include <exception>
#include <matrix2d.h>

class MINIMUMSHARED_EXPORT DeltaGridException : public std::exception
{
public:
    DeltaGridException(const std::string &msg = "");
    virtual const char* what() const noexcept;
private:
    std::string message;
};

class MINIMUMSHARED_EXPORT DeltaGrid2D
{
public:
    DeltaGrid2D();
    virtual ~DeltaGrid2D();

    auto initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void;
    auto cleanGrid() -> void;

    auto distributeGauss(const SpacePoint& sp, unsigned int sigmaXNum = 1, unsigned int sigmaYNum = 1) -> void;
    auto distributeSigle(const SpacePoint& sp) -> void;
    auto distributeRect4(const SpacePoint& sp) -> void;

    auto consentrateInPoint(const DoubleMatrix &m, int v=4) const -> double;
    auto consentrateInPoint(const DoubleMatrix &m, double &dx, double &dy) const -> double;

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

    auto hx() const -> double;
    auto hy() const -> double;

    double **nodes() { return m_nodes; }
    double **nodes() const { return m_nodes; }

    bool *_rows;
    bool *_cols;

private:

    unsigned int _N;
    unsigned int _M;
    double _hx;
    double _hy;
    SpacePoint _p;

    double **m_nodes;
    unsigned int _rx;
    unsigned int _ry;
    unsigned int _minX;
    unsigned int _maxX;
    unsigned int _minY;
    unsigned int _maxY;
};

class MINIMUMSHARED_EXPORT DeltaGrid1D
{
public:
    DeltaGrid1D();
    virtual ~DeltaGrid1D();

    auto initGrid(unsigned int N, double hx) -> void;
    auto cleanGrid() -> void;

    auto distributeGauss(const SpacePoint& sp, unsigned int sigmaXNum = 1) -> void;
    auto distributeSigle(const SpacePoint& sp) -> void;
    auto distributeRect4(const SpacePoint& sp) -> void;

    auto consentrateInPoint(const DoubleVector &m, int v) const -> double;
    auto consentrateInPoint(const DoubleVector &m, double &dx) const -> double;

    auto isCenter(const SpaceNodePDE &sn) const -> bool;
    auto isCenter(unsigned int n) const -> bool;

    auto isContains(const SpaceNodePDE &sn) const -> bool;
    auto isContains(unsigned int n) const -> bool;

    auto weight(const SpaceNodePDE &sn) const -> double;
    auto weight(unsigned int n) const -> double;

    auto rx() const -> unsigned int;

    auto p() const -> const SpacePoint&;
    auto p() -> SpacePoint&;

    auto minX() const -> unsigned int;
    auto maxX() const -> unsigned int;

    double *nodes() { return m_nodes; }
    double *nodes() const { return m_nodes; }

    auto hx() const -> double;

private:

    unsigned int _N;
    double _hx;
    SpacePoint _p;

    double *m_nodes;
    unsigned int _rx;
    unsigned int _minX;
    unsigned int _maxX;
};

#endif // DELTA_GRID_H
