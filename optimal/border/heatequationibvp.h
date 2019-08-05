#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "border_global.h"

class BORDERSHARED_EXPORT HeatEquationIBVP : public IHeatEquationIBVP
{
public:
    static void Main(int argc, char *argv[]);

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleVector&, const TimeNodePDE&) const;
    virtual void layerInfo(const DoubleMatrix&, const TimeNodePDE&) const;

    double U(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

protected:
    double lambda() const { return 0.5; }
};

class BORDERSHARED_EXPORT FinalHeatEquationIBVP : public IFinalHeatEquationIBVP
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
