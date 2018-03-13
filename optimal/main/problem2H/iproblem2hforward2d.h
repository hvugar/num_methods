#ifndef IPROBLEM2HFORWARD2D_H
#define IPROBLEM2HFORWARD2D_H

#include <grid/hibvp.h>
#include <vector>

using namespace std;

struct Parameter
{
    unsigned int No;
    unsigned int Nc;
    unsigned int Ns;
    DoubleMatrix k;
    DoubleMatrix z;
    std::vector<SpacePoint> xi;
    std::vector<SpacePoint> eta;
    std::vector<SpacePoint> theta;
    DoubleVector q;
};

struct ExtendedSpacePointNode
{
    SpacePoint pt;
    unsigned int id;
    unsigned int i;
    unsigned int j;
    double x;
    double y;
    double w;
};

struct ObservationPointNode : public ExtendedSpacePointNode {};

struct ObservationDeltaNode : public ExtendedSpacePointNode {};

struct ControlPointNode : public ExtendedSpacePointNode {};

struct ControlDeltaNode : public ExtendedSpacePointNode {};

class IProblem2HForward2D : public HyperbolicIBVP
{
public:
    void calculateMVD(DoubleMatrix &u) const;
    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;

    Parameter mParameter;

protected:
    void extendObservationPoint(const SpacePoint &xi, std::vector<ObservationPointNode> &ons, unsigned int j) const;
    void extendContrlDeltaPoint(const SpacePoint &eta, std::vector<ControlDeltaNode> &cps, unsigned int id) const;
    void distributeDelta(const SpacePoint &pt, std::vector<ExtendedSpacePointNode> &nodes, unsigned int id) const;

protected:
    virtual double initial1(const SpaceNodePDE &sn) const;
    virtual double initial2(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
};

#endif // IPROBLEM2HFORWARD2D_H
