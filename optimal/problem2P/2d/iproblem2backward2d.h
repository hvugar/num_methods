#ifndef IPROBLEM2BACKWARD2D_H
#define IPROBLEM2BACKWARD2D_H

#include "iproblem2pibvp2d.h"

using namespace std;

class IProblem2Backward2D : public IProblem22DPIBVP
{
public:
    virtual ~IProblem2Backward2D();

    void calculateMVD(DoubleMatrix &p, vector<ExtendedSpaceNode2D> &info, bool use = true);
    void calculateMVD1(DoubleMatrix &p, vector<ExtendedSpaceNode2D> &info, bool use = true);
    virtual void layerInfo(const DoubleMatrix &p, unsigned int layerNumber) const = 0;

    virtual double penalty(unsigned int i, const TimeNodePDE &tn) const = 0;

    void setPenaltyCoefficient(double r);
    double penaltyCoefficient() const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

protected:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int j, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    bool checkDelta(double delta) const;

    void extendControlPoint0(const SpaceNodePDE &cp, std::vector<ControlNode> &ops, unsigned int i) const;
    void extendControlPoint1(const SpaceNodePDE &cp, std::vector<ControlNode> &ops, unsigned int i) const;
    void extendControlPoint2(const SpaceNodePDE &cp, std::vector<ControlNode> &ops, unsigned int i) const;
    void extendControlPoint3(const SpaceNodePDE &cp, std::vector<ControlNode> &ops, unsigned int i) const;
    void extendControlPoint4(const SpaceNodePDE &cp, std::vector<ControlNode> &ops, unsigned int i) const;

    void extendObservationDeltaPoint0(const SpaceNodePDE &op, std::vector<ObservationDeltaNode> &ops, unsigned int j) const;
    void extendObservationDeltaPoint1(const SpaceNodePDE &op, std::vector<ObservationDeltaNode> &ons, unsigned int j) const;
    void extendObservationDeltaPoint2(const SpaceNodePDE &op, std::vector<ObservationDeltaNode> &ons, unsigned int j) const;
    void extendObservationDeltaPoint3(const SpaceNodePDE &op, std::vector<ObservationDeltaNode> &ons, unsigned int j) const;
    void extendObservationDeltaPoint4(const SpaceNodePDE &op, std::vector<ObservationDeltaNode> &ons, unsigned int j) const;

protected:
    double r;

#ifdef USE_ADDITIONAL_FUNCTIONS
protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
#endif
};

#endif // IPROBLEM2BACKWARD2D_H
