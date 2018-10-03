#ifndef PROBLEM22DEX2_H
#define PROBLEM22DEX2_H

#include "../abs/abstractproblem22d.h"
#include <imaging.h>

class Problem2Forward2DEx2 : public IProblem2Forward2D
{
protected:
    virtual double initial(const SpaceNodePDE &) const { return 0.0; }
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual void layerInfo(const DoubleMatrix &u, unsigned int layerNumber) const;
};

class Problem2Backward2DEx2 : public IProblem2Backward2D
{
protected:
    virtual double initial(const SpaceNodePDE &) const { return 0.0; }
    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const { return NAN; }
    virtual double f(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    virtual double h(const SpaceNodePDE &) const { return 0.0; }
    virtual double g1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g2(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g3(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double g4(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
};

class Problem22DEx2 : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);
};

#endif // PROBLEM22DEX2_H
