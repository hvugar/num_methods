#ifndef PROBLEM1P_SOLVER_H
#define PROBLEM1P_SOLVER_H

#include "problem1p_global.h"
#include <function.h>
#include <gradient.h>
#include <printer.h>
#include <grid/pibvp.h>

namespace p1p
{

class ProblemSolver;

struct EquationParameters
{
    double thermalDiffusivity = 1.0;
    double thermalConductivity0 = 0.001;
    double thermalConductivity1 = 0.1;
    double thermalConductivity2 = 0.01;

    unsigned int L;
    double *k;
    double *z;
    SpacePoint *eta;

    double initialTemperature;
    double environmentTemperature;

    double *v;
};

class PROBLEM1PSHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;
};

class PROBLEM1PSHARED_EXPORT HeatEquationFBVP : public IHeatEquationFBVP
{
public:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;
};

class PROBLEM1PSHARED_EXPORT ProblemSolver : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char* argv[]);

    ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX);

    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    virtual double fx(const DoubleVector &x) const;

    void setTimeDimension(const Dimension &timeDimension);
    void setSpaceDimensionX(const Dimension &spaceDimensionX);

protected:
    HeatEquationIBVP forward;
    HeatEquationFBVP backward;

    EquationParameters params;
    Dimension _timeDimension;
    Dimension _spaceDimensionX;

    DoubleVector u;
    DoubleVector U;
    DoubleVector p0;

    ProblemSolver *const_this;

    friend class HeatEquationIBVP;
    friend class HeatEquationFBVP;
};

};

#endif // PROBLEM1P_SOLVER_H
