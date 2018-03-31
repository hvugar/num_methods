#ifndef CPROBLEM2BACKWARD2D_H
#define CPROBLEM2BACKWARD2D_H

#include "iproblem2backward2d.h"

/*
 * p(x,y,t)=x*x+y*y+t;
 */
class CProblem2Backward2D : public IProblem2Backward2D
{
public:
    static void Main(int argc, char* argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual double penalty(unsigned int i UNUSED_PARAM, const TimeNodePDE &tn UNUSED_PARAM) const { return 10.0; }

protected:
    virtual double g1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g2(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g3(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double g4(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    virtual double P(double x, double y, double t) const;
};

#endif // CPROBLEM2BACKWARD2D_H
