#ifndef HYPERBOLICIBVP_H
#define HYPERBOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The IHyperbolicIBVP class
 * u_tt(x,t) = a^2u_xx(x,t) + f(x,t), t in (0,T], x in (0,l)
 */
class MINIMUMSHARED_EXPORT IHyperbolicIBVP : public InitialBoundaryValueProblemPDE
{
public:
    virtual ~IHyperbolicIBVP();
protected:
    virtual double initial1(const SpaceNodePDE &sn) const = 0;
    virtual double initial2(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    
public:
    void calculateU1(DoubleVector &u, double a=1.0, double lambda=0.25) const;
    void calculateU2(DoubleVector &u, double a=1.0) const;
    void calculateU3(DoubleMatrix &u, double a=1.0, double lambda=0.25) const;
};

/**
 * @brief The HyperbolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t), t in (0,T], x in (0, l),
 * u(x,0)   = fi1(x), x in [0,l],
 * u_t(x,0) = fi2(x), x in [0,l],
 * u(0,t)   = m1(t),  t in (0,T],
 * u(1,t)   = m2(t),  t in (0,T].
 */

class MINIMUMSHARED_EXPORT HyperbolicIBVP : protected IHyperbolicIBVP
{
protected:
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);

    void gridMethod0(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod1(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod2(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
};

#endif // HYPERBOLICIBVP_H
