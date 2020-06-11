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
    DeltaGrid2D(const DeltaGrid2D &dg);
    virtual ~DeltaGrid2D();

    DeltaGrid2D& operator =(const DeltaGrid2D &dg);

    auto initGrid(unsigned int N, double hx, unsigned int M, double hy) -> void;
    auto initGrid(const Dimension &dimensionX, const Dimension &dimensionY) -> void;
    auto cleanGrid() -> void;
    auto reset() -> void;
    auto resetAll() -> void;

    auto gaussWeight(const SpacePoint &sp, const SpacePoint &mu, double sigmaX, double sigmaY) const -> double;
    auto distributeGauss(const SpacePoint &sp, unsigned int nodeX_per_sigmaX = 1, unsigned int nodeY_per_sigmaY = 1) -> void;
    auto distributeGauss1(const SpacePoint &sp, unsigned int nodeX_per_sigmaX = 1, unsigned int nodeY_per_sigmaY = 1) -> void;
    auto lumpPointGauss(const DoubleMatrix &m) const -> double;
    auto lumpPointGauss(const DoubleMatrix &u, double &dx, double &dy) const -> double;
    auto lumpPointGauss1(const DoubleMatrix &u, double &dx, double &dy) const -> double;

    auto distributeSigle(const SpacePoint& sp) -> void;
    auto distributeRect4(const SpacePoint& sp) -> void;

    auto consentrateInPoint(const DoubleMatrix &u, unsigned int v=4) const -> double;
    auto derivativesInPoint(const DoubleMatrix &u, double &dx, double &dy, unsigned int v=4) const -> void;
    auto consentrateInPoint(const DoubleMatrix &u, double &dx, double &dy, unsigned int v=4) const -> double;

    auto isCenter(const SpaceNodePDE &sn) const -> bool;
    auto isCenter(unsigned int n, unsigned int m) const -> bool;

    auto isContains(const SpaceNodePDE &sn) const -> bool;
    auto isContains(unsigned int n, unsigned int m) const -> bool;

    auto weight(const SpaceNodePDE &sn) const -> double;
    auto weight(unsigned int n, unsigned int m) const -> double;

    auto rx() const -> unsigned int;
    auto ry() const -> unsigned int;

    auto p() const -> const SpacePoint &;
    auto p() -> SpacePoint &;

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

    double onFlyWeight(const SpacePoint &sp, const SpacePoint &mu, size_t n) const;

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
    DeltaGrid1D(const DeltaGrid1D &dg);
    virtual ~DeltaGrid1D();

    DeltaGrid1D& operator =(const DeltaGrid1D &dg);

    auto initGrid(unsigned int N, double hx) -> void;
    auto initGrid(const Dimension &dimension) -> void;
    auto cleanGrid() -> void;
    auto reset() -> void;
    auto resetAll() -> void;

    auto distributeGauss(const SpacePoint &sp, unsigned int nodeX_per_sigmaX = 1) -> void;
    auto lumpPointGauss(const DoubleVector &u) const -> double;
    auto lumpPointGauss(const DoubleVector &u, double &dx) const -> double;
    auto lumpPointGauss1(const DoubleVector &u, double &dx) const -> double;

    auto distributeSigle(const SpacePoint &sp) -> void;
    auto distributeRect4(const SpacePoint &sp) -> void;

    auto consentrateInPoint(const DoubleVector &u, unsigned int v=4) const -> double;
    auto derivativesInPoint(const DoubleVector &u, double &dx, unsigned int v=4) const -> void;
    auto consentrateInPoint(const DoubleVector &u, double &dx, unsigned int v=4) const -> double;

    auto isCenter(const SpaceNodePDE &sn) const -> bool;
    auto isCenter(unsigned int n) const -> bool;

    auto isContains(const SpaceNodePDE &sn) const -> bool;
    auto isContains(unsigned int n) const -> bool;

    auto weight(const SpaceNodePDE &sn) const -> double;
    auto weight(unsigned int n) const -> double;

    auto rx() const -> unsigned int;

    auto p() const -> const SpacePoint &;
    auto p() -> SpacePoint &;

    auto minX() const -> unsigned int;
    auto maxX() const -> unsigned int;

    double *nodes() { return m_nodes; }
    double *nodes() const { return m_nodes; }

    double *der_x() const { return m_der_x; }

    auto hx() const -> double;

    double onFlyWeight(const SpaceNodePDE &sn, const SpacePoint &mu, size_t n) const;

private:
    unsigned int _N;
    double _hx;
    SpacePoint _p;

    double *m_nodes;
    double *m_der_x;
    unsigned int _rx;
    unsigned int _minX;
    unsigned int _maxX;
};

class MINIMUMSHARED_EXPORT DeltaFunction
{
public:
    static double gaussian(double p, double m, double sigma);
    static double gaussian(const SpacePoint &p, const SpacePoint &m, const SpacePoint &sigma);
    static double gaussian(double p, double m, double sigma, size_t k);
    static double gaussian(const SpacePoint &p, const SpacePoint &m, const SpacePoint &sigma, size_t k);
    static double sinusoid(double p);
    static double sinusoid(const SpacePoint &p);

    static double lumpedPoint2(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY);
    static double lumpedPoint2(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny);
    static double lumpedPoint3(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY);
    static double lumpedPoint3(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny);
    static double lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY);
    static double lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, double hx, double hy, unsigned int Nx, unsigned int Ny);
    static double lumpedPoint4(const DoubleMatrix &m, const SpacePoint &p, const Dimension &dimensionX, const Dimension &dimensionY,
                               SpacePoint &d);
    static double lumpedPointG(const DoubleMatrix &u, const SpacePoint &m, const Dimension &dimensionX, const Dimension &dimensionY, unsigned int nps, size_t k);
};

#endif // DELTA_GRID_H
