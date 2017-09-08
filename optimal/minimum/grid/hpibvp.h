#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "pibvp.h"

class HeatEquationIBVP : protected InitialBoundaryValueProblemPDE
{
public:
    virtual ~HeatEquationIBVP();
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, double a = 1.0, SweepMethodDirection direction = ForwardSweep);
};

#endif // HEATEQUATIONIBVP_H
