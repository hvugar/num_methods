#ifndef PROBLEM2FORWARD2DEX4_H
#define PROBLEM2FORWARD2DEX4_H

#include "../iproblem2forward2d.h"

class Problem2Forward2DEx4 : public IProblem2Forward2D
{
public:
    virtual ~Problem2Forward2DEx4() {}
    double fi;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;
};

#endif // PROBLEM2FORWARD2DEX4_H
