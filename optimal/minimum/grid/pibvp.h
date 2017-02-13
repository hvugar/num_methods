#ifndef PARABOLICIBVP_H
#define PARABOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The ParabolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t),
 *          t in (0,T], x in (0, l),
 * u(x,0) = fi(x), x in [0,l],
 * u(0,t) = m1(t), t in (0,T],
 * u(l,t) = m2(t), t in (0,T].
 */
class MINIMUMSHARED_EXPORT ParabolicIBVP : protected InitialBoundaryValueProblemPDE
{
public:
    virtual ~ParabolicIBVP();

protected:
    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);

    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u);
    //void calculateN2R2LD(DoubleMatrix &u);
    void calculateN4L2RD(DoubleMatrix &u);
    //void calculateN4R2LD(DoubleMatrix &u);
    //void calculateN6L2RD(DoubleMatrix &u);
    //void calculateN6R2LD(DoubleMatrix &u);

    void calculateMVD(DoubleMatrix &u) const;
};

#endif // PARABOLICIBVP_H
