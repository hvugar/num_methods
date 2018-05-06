#ifndef IPROBLEM2HBACKWARD2D_H
#define IPROBLEM2HBACKWARD2D_H

#include "iproblem2h2d.h"
#include "iproblem2h2d_ifunctional.h"

namespace IProblem2H
{

class IProblem2HBackward2D : public IHyperbolicIBVP
{
public:
    void calculateMVD(DoubleMatrix &p, std::vector<SpacePointInfo> &p_info, bool use, const std::vector<SpacePointInfo> &u_info) const;
private:
    void calculateMVD_D(DoubleMatrix &p, std::vector<SpacePointInfo> &p_info, bool use, const std::vector<SpacePointInfo> &u_info) const;
    void calculateMVD_N(DoubleMatrix &p, std::vector<SpacePointInfo> &p_info, bool use, const std::vector<SpacePointInfo> &u_info) const;

public:
    virtual void layerInfo(const DoubleMatrix &p, unsigned int layerNumber) const;

    void add2Info(const DoubleMatrix &p, std::vector<SpacePointInfo> &info, const std::vector<ExtendedSpacePointNode> &cntPointNodes, unsigned int ln, double hx, double hy, unsigned int rows, unsigned int cols) const;

    EquationParameter mEquParameter;
    OptimizeParameter mOptParameter;

    IFunctional *ifunc;
    DoubleMatrix UT;
    DoubleMatrix UTt;
protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

}

#endif // IPROBLEM2HBACKWARD2D_H
