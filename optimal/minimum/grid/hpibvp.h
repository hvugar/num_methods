#ifndef HEATEQUATIONIBVP_H
#define HEATEQUATIONIBVP_H

#include "pibvp.h"

class MINIMUMSHARED_EXPORT HeatEquationIBVP : public InitialBoundaryValueProblemPDE
{
public:
    virtual ~HeatEquationIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryType boundary = Unused) const = 0;
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

#endif // HEATEQUATIONIBVP_H
