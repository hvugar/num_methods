#ifndef IPROBLEM2FORWARD2D_H
#define IPROBLEM2FORWARD2D_H

#define USE_OTHER_FUNCTIONS_F

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include "problem2setting.h"
#include <time.h>

using namespace std;

struct ObservationNode
{
    SpaceNodePDE xi;
    unsigned int j;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};

struct ControlDeltaNode
{
    SpaceNodePDE eta;
    unsigned int i;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};

class IProblem2Forward2D : public IParabolicIBVP
{
public:
    void setSettings(P2Setting s);

    void calculateMVD(std::vector<DoubleMatrix> &u) const;
    void calculateMVD(DoubleMatrix &u) const;

    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    bool checkDelta(double delta) const;

    void extendContrlDeltaPoint(const SpaceNodePDE cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;
    void extendContrlDeltaPoint1(const SpaceNodePDE cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;
    void extendContrlDeltaPoint2(const SpaceNodePDE cp, std::vector<ControlDeltaNode> &cps, unsigned int i) const;

    void extendObservationPoint1(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;

    void extendObservationPoint(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;

protected:
    P2Setting setting;
};

#endif // IPROBLEM2FORWARD2D_H
