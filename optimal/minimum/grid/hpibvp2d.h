#ifndef HPIBVP2D_H
#define HPIBVP2D_H

#include "pibvp.h"

class HeatEquationIBVP2D : public InitialBoundaryValueProblemPDE
{
public:
    virtual ~HeatEquationIBVP2D();

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double env0(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double env1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    /**
     * @brief calculateU
     * u_t = a^2(u_xx + u_yy) - alpha(u(x,y,t)-theta(x,y,t)) + f(x,y,t);
     * u(x,y,0) = \phi(x,y);
     * u_n(x,y,t) = lambda(u(x,y,t)-m(x,y,t));
     * @param u
     * @param a
     * @param alpha
     * @param lambda
     */
    void calculateU(DoubleMatrix &u, double a, double alpha, double lambda);
};

#endif // HPIBVP2D_H
