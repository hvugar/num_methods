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

public:
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);

    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u);
    void calculateN2R2LD(DoubleMatrix &u);

    void calculateN4L2RD(DoubleMatrix &u);
    void calculateN4R2LD(DoubleMatrix &u);

    void calculateN6L2RD(DoubleMatrix &u);
    void calculateN6R2LD(DoubleMatrix &u);
};

#endif // PARABOLICIBVP_H
