#ifndef CPROBLEM2FORWARD2D_H
#define CPROBLEM2FORWARD2D_H

#include "iproblem2forward2d.h"

class CProblem2Forward2D : public IProblem2Forward2D
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}
public:
    double U(double x, double y, double t) const;
};

#endif // CPROBLEM2FORWARD2D_H
