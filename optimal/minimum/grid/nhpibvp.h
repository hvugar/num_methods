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
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const = 0;

    virtual double theta0(const TimeNode &tn) const = 0;
    virtual double theta1(const TimeNode &tn) const = 0;
    virtual double theta2(const TimeNode &tn) const = 0;

public:
    void calculateGM1(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void calculateGM2(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
};

#endif // NEWTONHEATEQUATION_H
