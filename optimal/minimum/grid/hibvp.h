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
};

class MINIMUMSHARED_EXPORT CCIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CCIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void calculateD1V1(DoubleVector &u, double a) const;
    void calculateD1V2(DoubleVector &u, double a, double lambda=0.25) const;
    void calculateD2V1(DoubleMatrix &u, double a) const;
    void calculateD2V2(DoubleMatrix &u, double a, double lambda=0.25) const;
};

class MINIMUMSHARED_EXPORT CC1IHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CC1IHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicitCalculateD2V1(DoubleMatrix &m, double a, double alpha) const;
    void implicitCalculateD2V1(DoubleMatrix &m, double a, double alpha) const;

    void calculateD2V1(DoubleMatrix &u, double a, double sigma) const;
    void calculateD2V2(DoubleMatrix &u, double a, double lambda=0.25) const;

private:
    auto calculateInitialD2(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, unsigned int M, double hx, double hy, double ht, double a, double sigma) const -> void;
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
public:
    virtual ~HyperbolicIBVP();

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
