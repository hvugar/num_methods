#ifndef NEWTONHEATEQUATION_H
#define NEWTONHEATEQUATION_H

#include "pibvp.h"

class MINIMUMSHARED_EXPORT NewtonHeatEquation : public IParabolicIBVP
{
public:
    double lambda0;
    double lambda1;
    double lambda2;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double theta0(const TimeNodePDE &tn) const = 0;
    virtual double theta1(const TimeNodePDE &tn) const = 0;
    virtual double theta2(const TimeNodePDE &tn) const = 0;

public:
    void calculateGM1(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void calculateGM2(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void calculateGM3(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
};

#endif // NEWTONHEATEQUATION_H
