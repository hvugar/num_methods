#ifndef IPROBLEM2HFORWARD2D_H
#define IPROBLEM2HFORWARD2D_H

#include "iproblem2h2d.h"

namespace IProblem2H
{

class IProblem2HForward2D : public IHyperbolicIBVP
{
public:
    void calculateMVD(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
private:
    void calculateMVD_N(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
    void calculateMVD_D(DoubleMatrix &u, DoubleMatrix &ut, vector<ExtendedSpaceNode2DH> &info, bool use) const;
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
    void add2Info(const DoubleMatrix &u, vector<ExtendedSpaceNode2DH> &info, unsigned int ln) const;
};

}

#endif // IPROBLEM2HFORWARD2D_H
