#ifndef PROBLEM22DEX3_H
#define PROBLEM22DEX3_H

#include "../abs/abstractproblem22d.h"
#include <imaging.h>

class Problem2Forward2DEx3;
class Problem2Backward2DEx3;
class Problem22DEx3;

class Problem2Forward2DEx3 : public IProblem2Forward2D
{
protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}
};

class Problem2Backward2DEx3 : public IProblem2Backward2D
{
public:
    Problem22DEx3 *p22dEx3;

protected:
    virtual double initial(const SpaceNodePDE &) const;
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double h(const SpaceNodePDE &) const { return 0.0; }
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
};

class Problem22DEx3 : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM22DEX3_H
