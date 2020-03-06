#ifndef PROBLEM1P_SOLVER_H
#define PROBLEM1P_SOLVER_H

#define SIGMA 0.5
#define EXAMPLE_LEFT_BORDER_ROBIN
//#define EXAMPLE_LEFT_BORDER_DIRICHLET
//#define EXAMPLE_FXT

#include "problem1p_global.h"
#include <function.h>
#include <gradient.h>
#include <printer.h>
#include <grid/pibvp.h>

#include <iostream>
#include <functional>
#include <memory>

using namespace std;
using namespace std::placeholders;

namespace p1p
{

class ProblemSolver;

struct EquationParameters
{



    unsigned int L;
    double *k;
    double *z;
    SpacePoint *eta;


};

class PROBLEM1PSHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;

protected:
    virtual double weight() const { return SIGMA; }

public:
    virtual Dimension timeDimension() const;
    virtual Dimension spaceDimensionX() const;
    virtual Dimension spaceDimensionY() const;
    virtual Dimension spaceDimensionZ() const;
};

class PROBLEM1PSHARED_EXPORT HeatEquationFBVP : public IHeatEquationFBVP
{
public:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    ProblemSolver *solver;

protected:
    virtual double weight() const { return SIGMA; }

public:
    virtual Dimension timeDimension() const;
    virtual Dimension spaceDimensionX() const;
    virtual Dimension spaceDimensionY() const;
    virtual Dimension spaceDimensionZ() const;
};

class PROBLEM1PSHARED_EXPORT ProblemSolver : public RnFunction, public IGradient
{
public:
    static void Main(int argc, char* argv[]);

    ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX);

    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    virtual double fx(const DoubleVector &x) const;

    const Dimension &timeDimension() const { return _timeDimension; }
    void setTimeDimension(const Dimension &timeDimension);

    const Dimension &spaceDimensionX() const { return _spaceDimensionX; }
    void setSpaceDimensionX(const Dimension &spaceDimensionX);

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleVector &u, const TimeNodePDE &tn) const;

    double integral(const DoubleVector &x) const;

protected:
    HeatEquationIBVP forward;
    HeatEquationFBVP backward;

    double environmentTemperature = 0.0;
    double initialTemperature = 0.5;

    double thermalDiffusivity = 1.0;
    double thermalConvection = 0.001;
    double thermalConductivity1 = 0.1;
    double thermalConductivity2 = 0.01;

    double *heat_power = nullptr;

    EquationParameters params;

    DoubleVector U;
    DoubleVector V;
    DoubleVector p0;
    DoubleVector p1;
    DoubleVector p2;
    DoubleVector p0x;

    ProblemSolver *const_this;

private:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;

    friend class HeatEquationIBVP;
    friend class HeatEquationFBVP;
};

};

#endif // PROBLEM1P_SOLVER_H
