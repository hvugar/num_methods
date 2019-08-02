#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "border_global.h"

class HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    static void Main(int argc, char *argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;
};

class FinalHeatEquationIBVP : public IFinalHeatEquationIBVP
{
public:
    static void Main(int argc, char *argv[]);

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;
};

#endif // HEATEQUATIONIBVP_H
