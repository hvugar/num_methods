#ifndef HEATEQAUATIONIBVP_H
#define HEATEQAUATIONIBVP_H

#include "pibvp.h"

class HeatEqauationIBVP : protected InitialBoundaryValueProblemPDE
{
public:
    virtual ~HeatEqauationIBVP();
protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, double a = 1.0, SweepMethodDirection direction = ForwardSweep);
};

#endif // HEATEQAUATIONIBVP_H
