#ifndef DELTA_GRID_H
#define DELTA_GRID_H

#include "grid/grid.h"
#include <cmath>
#include <exception>
#include <matrix2d.h>
#include "exceptions.h"

class MINIMUMSHARED_EXPORT DeltaGrid2D
{
public:
    DeltaGrid2D();
    DeltaGrid2D(const DeltaGrid2D&);
    virtual ~DeltaGrid2D();

    DeltaGrid2D& operator =(const DeltaGrid2D& vector);


    auto initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void;
    auto initGrid(const Dimension &dimensionX, const Dimension &dimensionY) -> void;
    auto cleanGrid() -> void;
    auto reset() -> void;
    auto resetAll() -> void;

    auto distributeGauss(const SpacePoint& sp, unsigned int nodeX_per_sigmaX = 1, unsigned int nodeY_per_sigmaY = 1) -> void;
    auto lumpPointGauss(const DoubleMatrix &m) const -> double;
    auto lumpPointGauss(const DoubleMatrix &m, double &mx, double &my) const -> double;

    auto distributeSigle(const SpacePoint& sp) -> void;
    auto distributeRect4(const SpacePoint& sp) -> void;

    auto consentrateInPoint(const DoubleMatrix &m, unsigned int v=4) const -> double;
    auto derivativesInPoint(const DoubleMatrix &m, double &dx, double &dy, unsigned int v=4) const -> void;
    auto consentrateInPoint(const DoubleMatrix &m, double &dx, double &dy, unsigned int v=4) const -> double;

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

    double **der_x() const { return m_der_x; }
    double **der_y() const { return m_der_y; }

private:

    unsigned int _N;
    unsigned int _M;
    double _hx;
    double _hy;
    SpacePoint _p;

    double **m_nodes;
    double **m_der_x;
    double **m_der_y;
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
