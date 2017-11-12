#ifndef IPROBLEM2FORWARD2D_H
#define IPROBLEM2FORWARD2D_H

#include <grid/pibvp.h>
#include <vector>
#include <printer.h>

using namespace std;

class IProblem2Forward2D : public IParabolicIBVP
{
public:
    IProblem2Forward2D();
    virtual ~IProblem2Forward2D() {}
    void setSettings(double a, double lambda0, double lambda, double theta, unsigned int Lc, unsigned int Lo);
    void calculateMVD(DoubleMatrix &u) const;

    void calculateMVD1(DoubleMatrix &u) const;
    void calculateMVD1X(DoubleMatrix &u, DoubleMatrix &uh, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn, double ht, unsigned int *dmx) const;
    void calculateMVD1Y(DoubleMatrix &u, DoubleMatrix &uh, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn, double ht, unsigned int *dmy) const;

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

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;

public:
    double U(double x, double y, double t) const;

private:
    double a;
    double lambda0;
    double lambda;
    double theta;

    unsigned int Lc;
    unsigned int Lo;

    DoubleMatrix k;
    DoubleMatrix z;
    vector<SpaceNodePDE> xi;
    vector<SpaceNodePDE> eta;
};

#endif // IPROBLEM2FORWARD2D_H
