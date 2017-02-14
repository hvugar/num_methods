#ifndef BACKWARDPARABOLICIBVP_H
#define BACKWARDPARABOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The BackwardParabolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t),
 *          t in [0,T), x in (0, l),
 * u(x,T) = fi(x), x in [0,l],
 * u(0,t) = m1(t), t in [0,T),
 * u(l,t) = m2(t), t in [0,T).
 */
class MINIMUMSHARED_EXPORT BackwardParabolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:

    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;
    void calculate(DoubleCube &u, SweepMethodDirection direction = ForwardSweep) const;
};

#endif // BACKWARDPARABOLICIBVP_H
