#ifndef IPROBLEM2HFORWARD2D_H
#define IPROBLEM2HFORWARD2D_H

#include "iproblem2h2d.h"

namespace IProblem2H
{

class IProblem2HForward2D : public IHyperbolicIBVP
{
public:
    void calculateMVD(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
    void extendPoints(std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                      std::vector<ExtendedSpacePointNode> &qPointNodes, const Dimension &dimX, const Dimension &dimY) const;
    void findRowsCols(std::vector<unsigned int> &rows0, std::vector<unsigned int> &rows1, std::vector<unsigned int> &rows2,
                      std::vector<unsigned int> &cols0, std::vector<unsigned int> &cols1, std::vector<unsigned int> &cols2,
                      std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                      unsigned int N, unsigned int M) const;
    void initialLayers(DoubleMatrix &u00, DoubleMatrix &u05, DoubleMatrix &u10, vector<ExtendedSpaceNode2DH> &info, bool use,
                       std::vector<ExtendedSpacePointNode> &obsPointNodes, std::vector<ExtendedSpacePointNode> &cntDeltaNodes,
                       std::vector<ExtendedSpacePointNode> &qPointNodes,unsigned int N, unsigned int M,
                       double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda) const;
    void mallocVectors(double *&a1X, double *&b1X, double *&c1X, double *&d1X, double *&x1X, unsigned int sizeX,
                       double *&a1Y, double *&b1Y, double *&c1Y, double *&d1Y, double *&x1Y, unsigned int sizeY) const;
    void mallocVectorsX(double *&a1X, double *&b1X, double *&c1X, double *&d1X, double *&x1X, unsigned int sizeX, double k0, double k1) const;
    void mallocVectorsY(double *&a1Y, double *&b1Y, double *&c1Y, double *&d1Y, double *&x1Y, unsigned int sizeY, double k0, double k1) const;
    void initBorders(unsigned int N, unsigned int M, double hx, double hy, double ht, unsigned int l,
                     DoubleMatrix &u15, DoubleMatrix &u) const;
    void freeVectors(double *a1X, double *b1X, double *c1X, double *d1X, double *x1X,
                     double *a1Y, double *b1Y, double *c1Y, double *d1Y, double *x1Y) const;
    void prepareInfo(unsigned int N, std::vector<SpacePoint> points, std::vector<ExtendedSpaceNode2DH> &info,
                     unsigned int L, const Dimension &dimX, const Dimension &dimY) const;

private:
    void calculateMVD_N(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;

    void calculateMVD_D(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
    void calculateMVD_D1(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
    void calculateMVD_D2(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
    void calculateMVD_D3(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
public:
    void calculateMVD_N(DoubleMatrix &u, DoubleMatrix &ut) const;

    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;
    virtual void layerInfo(const DoubleMatrix &u, const DoubleMatrix &ut, unsigned int layerNumber) const;

    EquationParameter mEquParameter;
    OptimizeParameter mOptParameter;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    void add2Info(const DoubleMatrix &u, vector<ExtendedSpaceNode2DH> &info, const std::vector<ExtendedSpacePointNode> &obsPointNodes,
                  unsigned int ln, double hx, double hy, unsigned int rows, unsigned int cols) const;
};

}

#endif // IPROBLEM2HFORWARD2D_H
