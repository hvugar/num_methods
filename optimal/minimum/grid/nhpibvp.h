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
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double theta0(const TimeNodePDE &tn) const = 0;
    virtual double theta1(const TimeNodePDE &tn) const = 0;
    virtual double theta2(const TimeNodePDE &tn) const = 0;

public:
    void calculateGM1(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void calculateGM2(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);
    void calculateGM3(DoubleVector &u, SweepMethodDirection direction = ForwardSweep);

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
};

/**
 * @brief The INewtonHeatEquation class
 * u_t(x,y,z,t) = a(x,y,z,t)div(grad u(x,y,z,t)) - lambda(x,y,z,t)[u(x,y,z,t)-theta(t)] + f(x,y,z,t),
 * u(x,y,z,0) = initial(x,y,z),
 * alpha(x,y,z,t)u(x,y,z,t)+beta(x,y,z,y)u/n(x,y,z,t)=border(x,y,z,t).
 */
class INewtonHeatEquation : public IParabolicIBVP
{
public:
    virtual ~INewtonHeatEquation() {}

protected:
    virtual auto a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
    virtual auto lambda(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
    virtual auto theta(const TimeNodePDE &tn) const -> double = 0;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
    virtual auto initial(const SpaceNodePDE &sn) const -> double = 0;
    virtual auto alpha(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
    virtual auto beta(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;
    virtual auto border(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double = 0;

    virtual auto layerInfo(DoubleVector &u, unsigned int ln) const -> void;
    virtual auto layerInfo(DoubleMatrix &u, unsigned int ln) const -> void;
public:
    auto calculate() const -> void;
};

#endif // NEWTONHEATEQUATION_H
