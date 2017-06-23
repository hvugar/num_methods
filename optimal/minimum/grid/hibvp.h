#ifndef HYPERBOLICIBVP_H
#define HYPERBOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The HyperbolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t)
 *          t in (0,T], x in (0, l),
 * u(x,0)   = fi1(x), x in [0,l],
 * u_t(x,0) = fi2(x), x in [0,l],
 * u(0,t) = m1(t), t in (0,T],
 * u(1,t) = m2(t), t in (0,T].
 */

class MINIMUMSHARED_EXPORT HyperbolicIBVP : protected InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial1(const SpaceNode &sn) const = 0;
    virtual double initial2(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const = 0;

public:
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);

    void gridMethod0(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod1(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod2(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
};

#endif // HYPERBOLICIBVP_H
