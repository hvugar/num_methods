#ifndef IPROBLEM2FORWARD2D_H
#define IPROBLEM2FORWARD2D_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>
#include <utils/random.h>
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

class IProblem2Forward2D : public IParabolicIBVP
{
public:
    IProblem2Forward2D();
    virtual ~IProblem2Forward2D() {}

    void setSettings(P2Setting s);

    void calculateMVD(std::vector<DoubleMatrix> &u) const;

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

public:
    double delta(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i, unsigned int source = 10) const;
    double delta1(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta2(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta3(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    double delta4(const SpaceNodePDE &sn, const SpaceNodePDE &eta, unsigned int i = 0) const;
    bool checkDelta(double delta) const;

    void extendObservationPoint1(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;
    void extendObservationPoint(const SpaceNodePDE op, std::vector<ObservationNode> &ops, unsigned int j) const;

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;

public:
    double U(double x, double y, double t) const;

private:
    //double a;
    //double lambda0;
    //double lambda;
    //double theta;

    //unsigned int Lc;
    //unsigned int Lo;

    //DoubleMatrix k;
    //DoubleMatrix z;
    //vector<SpaceNodePDE> xi;
    //vector<SpaceNodePDE> eta;

    P2Setting setting;
};

#endif // IPROBLEM2FORWARD2D_H
