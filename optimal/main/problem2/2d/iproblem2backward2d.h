#ifndef IPROBLEM2BACKWARD2D_H
#define IPROBLEM2BACKWARD2D_H

#include "iproblem2pibvp2d.h"

using namespace std;

struct ControlNode
{
    SpaceNodePDE eta;
    unsigned int i;
    unsigned int n;
    unsigned int m;
    double x;
    double y;
    double w;
};


class IProblem2Backward2D : public IProblem22DPIBVP
{
public:
    void setSettings(P2Setting s);

    void calculateMVD(std::vector<DoubleMatrix> &p);
    void calculateMVD(DoubleMatrix &p);

    virtual void layerInfo(const DoubleMatrix &p, unsigned int layerNumber) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int j, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    bool checkDelta(double delta) const;

    void extendControlPoint(const SpaceNodePDE cp, std::vector<ControlNode> &ops, unsigned int i) const;

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const = 0;
    virtual double h(const SpaceNodePDE &sn) const = 0;

protected:
    P2Setting setting;
};

#endif // IPROBLEM2BACKWARD2D_H
