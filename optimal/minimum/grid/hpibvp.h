#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "pibvp.h"

class MINIMUMSHARED_EXPORT HeatEquationIBVP : public InitialBoundaryValueProblemPDE
{
public:

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, double a = 1.0, SweepMethodDirection direction = ForwardSweep);

    /**
     * @brief
     * u_t = a^2 u_xx + f(x,t)
     * u(x,0) = 0
     * u_x(0,t) = u_x(l,t) = 0
     * @param u
     * @param a
     */
    void gridMethod1(DoubleVector &u, double a = 1.0);

};

class MINIMUMSHARED_EXPORT HeatEquationIBVP2D : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
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

#endif // HEATEQUATIONIBVP_H
