#ifndef PROBLEM2BACKWARD2D_H
#define PROBLEM2BACKWARD2D_H

#include "../iproblem2backward2d.h"

class AbstactProblem22D;

class Problem2Backward2D : public IProblem2Backward2D
{
public:
    virtual ~Problem2Backward2D() {}

    AbstactProblem22D *ap22d;
    DoubleMatrix *U;
    DoubleMatrix *u;
    vector<ExtendedSpaceNode2D> *info;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const;
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual double penalty(unsigned int i, const TimeNodePDE &tn) const;

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const;
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;
};

#endif // PROBLEM2BACKWARD2D_H
