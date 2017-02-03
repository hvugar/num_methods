#ifndef PARABOLICIBVP_H
#define PARABOLICIBVP_H

#include "ibvp.h"

/**
 * @brief The ParabolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t)
 * u(x,0) = fi(x),
 * u(0,t) = m1(t),
 * u(1,t) = m2(t).
 */
class MINIMUMSHARED_EXPORT ParabolicIBVP : protected InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(unsigned int n) const = 0;
    virtual double boundary(unsigned int m, BoundaryType boundary) const = 0;
    virtual double f(unsigned int n, unsigned int m) const = 0;
    virtual double a(unsigned int n, unsigned int m) const = 0;

    virtual double initial(unsigned int i, unsigned int j) const = 0;
    virtual double boundary(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double f(unsigned int i, unsigned int j, unsigned int k) const = 0;
    virtual double a(unsigned int i, unsigned int j, unsigned int k) const = 0;

    virtual double initial(const SpaceNode &sn) const = 0;
    virtual double boundary(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double f(const SpaceNode &sn, const TimeNode &tn) const = 0;
    virtual double a(const SpaceNode &sn, const TimeNode &tn) const = 0;

public:
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod2(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void calculateMVD(DoubleMatrix &u) const;
    void calculateMVD2(DoubleMatrix &u) const;

    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u);
    void calculateN2R2LD(DoubleMatrix &u);

    void calculateN4L2RD(DoubleMatrix &u);
    void calculateN4R2LD(DoubleMatrix &u);

    void calculateN6L2RD(DoubleMatrix &u);
    void calculateN6R2LD(DoubleMatrix &u);

    void calculateMVD1(DoubleMatrix &u) const;

    void calculateN2L2RD1(DoubleMatrix &u);
    void calculateN4L2RD1(DoubleMatrix &u);
};

#endif // PARABOLICIBVP_H
