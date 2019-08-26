#ifndef PROBLEM1P_SOLVER_H
#define PROBLEM1P_SOLVER_H

#include "problem1p_global.h"
#include <grid/pibvp.h>

namespace p1p
{

class PROBLEM1PSHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;
protected:
};

class PROBLEM1PSHARED_EXPORT ConjugateHeatEquationIBVP : public IFinalHeatEquationIBVP
{
public:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;
protected:
};

class PROBLEM1PSHARED_EXPORT Solver
{
protected:
    HeatEquationIBVP problem1;
    ConjugateHeatEquationIBVP problem2;

    unsigned int L;
};

};

#endif // PROBLEM1P_SOLVER_H
