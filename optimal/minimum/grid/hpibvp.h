#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "pibvp.h"

class HeatEquationIBVP : protected InitialBoundaryValueProblemPDE
{
public:
    virtual ~HeatEquationIBVP();
protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, double a = 1.0, SweepMethodDirection direction = ForwardSweep);
};

#endif // HEATEQUATIONIBVP_H
