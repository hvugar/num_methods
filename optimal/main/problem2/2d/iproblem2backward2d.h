#ifndef IPROBLEM2BACKWARD2D_H
#define IPROBLEM2BACKWARD2D_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include "problem2setting.h"

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


class IProblem2Backward2D : public IParabolicIBVP
{
public:
    IProblem2Backward2D();
    virtual ~IProblem2Backward2D() {}
    void setSettings(P2Setting s);

    void calculateMVD(std::vector<DoubleMatrix> &p);

    DoubleMatrix U;
    DoubleMatrix mu;
    DoubleMatrix uT;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &xi, unsigned int j, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;

    bool checkDelta(double delta) const;

    void extendControlPoint(const SpaceNodePDE cp, std::vector<ControlNode> &ops, unsigned int i) const;

    //double delta(const SpaceNodePDE &sn, unsigned int j, unsigned int source) const;

    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double h(const SpaceNodePDE &sn) const;

    virtual double P(double x, double y, double t) const;

private:
    P2Setting setting;
};

#endif // IPROBLEM2BACKWARD2D_H
