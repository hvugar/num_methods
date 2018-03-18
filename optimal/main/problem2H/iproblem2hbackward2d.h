#ifndef IPROBLEM2HBACKWARD2D_H
#define IPROBLEM2HBACKWARD2D_H

#include "iproblem2h2d.h"

class IProblem2HBackward2D : public IHyperbolicIBVP
{
public:
    void calculateMVD(DoubleMatrix &p) const;
    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;

    IProblem2H2D::Parameter mParameter;
    DoubleMatrix UT0;
    DoubleMatrix UT1;
protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    //void extendObservationDeltaPoint(const SpacePoint &xi, std::vector<IProblem2H2D::ObservationPointNode> &ons, unsigned int j) const;
    //void extendContrlPoint(const SpacePoint &eta, std::vector<IProblem2H2D::ControlDeltaNode> &cps, unsigned int id) const;
    //void distributeDelta(const SpacePoint &pt, std::vector<IProblem2H2D::ExtendedSpacePointNode> &nodes, unsigned int id) const;
};

#endif // IPROBLEM2HBACKWARD2D_H
