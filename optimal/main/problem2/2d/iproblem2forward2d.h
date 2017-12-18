#ifndef IPROBLEM2FORWARD2D_H
#define IPROBLEM2FORWARD2D_H

#include "iproblem2pibvp2d.h"

using namespace std;

class IProblem2Forward2D : public IProblem22DPIBVP
{
public:
    virtual ~IProblem2Forward2D();

    void calculateMVD(DoubleMatrix &u, vector<ExtendedSpaceNode2D> &info, bool use = true) const;
    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const = 0;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;

protected:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    bool checkDelta(double delta) const;

    void extendContrlDeltaPoint(const SpaceNodePDE &cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;
    void extendContrlDeltaPoint1(const SpaceNodePDE &cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;
    void extendContrlDeltaPoint2(const SpaceNodePDE &cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;

    void extendObservationPoint1(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;

    void extendObservationPoint(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;
    void extendObservationPoint(const SpaceNodePDE op, ExtendedSpaceNode2D &pi, unsigned int j) const;
    void extendObservationPoint(const SpaceNodePDE op, ExtendedSpaceNode2D &pi) const;
};

#endif // IPROBLEM2FORWARD2D_H
