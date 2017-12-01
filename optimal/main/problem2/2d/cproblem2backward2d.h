#ifndef CPROBLEM2BACKWARD2D_H
#define CPROBLEM2BACKWARD2D_H

#include "iproblem2backward2d.h"

class CProblem2Backward2D : public IProblem2Backward2D
{
public:
    DoubleMatrix U;
    DoubleMatrix mu;
    DoubleMatrix uT;
protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn UNUSED_PARAM) const;
    virtual double h(const SpaceNodePDE &sn) const;

    virtual double P(double x, double y, double t) const;
};

#endif // CPROBLEM2BACKWARD2D_H
