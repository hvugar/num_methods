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

class MINIMUMSHARED_EXPORT CcIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CcIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &u, double a) const;
    void implicit_calculate_D1V1(DoubleVector &u, double a) const;
    void implicit_calculate_D1V2(DoubleVector &u, double a, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &u, double a) const;
    void implicit_calculate_D2V1(DoubleMatrix &u, double a) const;
    void implicit_calculate_D2V2(DoubleMatrix &u, double a, double lambda=0.25) const;

    void implicit_calculate_D2V3(DoubleMatrix &u, double a) const;
    void implicit_calculate_D2V4(DoubleMatrix &u, double a) const;

private:
    void initial_calculate(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a) const;

    void initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a) const;
    void border1_calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const;
    void border2_calculate(DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn20) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT CdIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~CdIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &u, double a, double alpha) const;
    void implicit_calculate_D1V1(DoubleVector &u, double a, double alpha) const;
    void implicit_calculate_D1V2(DoubleVector &u, double a, double alpha, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const;
    void implicit_calculate_D2V1(DoubleMatrix &u, double a, double alpha) const;
    void implicit_calculate_D2V2(DoubleMatrix &u, double a, double alpha, double lambda=0.25) const;
    void implicit_calculate_D2V3(DoubleMatrix &u, double a, double alpha) const;

private:
    void initial_calculate(DoubleVector &u00, DoubleVector &u10, unsigned int N, double hx, double ht, double a, double sigma) const;

    void initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void border1_calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const;
    void border2_calculate(DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn20) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MINIMUMSHARED_EXPORT ConjugateCdIHyperbolicIBVP : public IHyperbolicIBVP
{
public:
    virtual ~ConjugateCdIHyperbolicIBVP();

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:
    void explicit_calculate_D1V1(DoubleVector &p, double a, double alpha) const;
    void implicit_calculate_D1V1(DoubleVector &p, double a, double alpha) const;
    void implicit_calculate_D1V2(DoubleVector &p, double a, double alpha, double lambda=0.25) const;

    void explicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const;
    void implicit_calculate_D2V1(DoubleMatrix &p, double a, double alpha) const;
    void implicit_calculate_D2V3(DoubleMatrix &p, double a, double alpha) const;
    void implicit_calculate_D2V2(DoubleMatrix &p, double a, double alpha, double lambda=0.25) const;

private:
    void initial_calculate(DoubleMatrix &u00, DoubleMatrix &u10, unsigned int N, double hx, unsigned int M, double hy, double ht, double a, double alpha) const;
    void border1_calculate(DoubleMatrix &u15, DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn15, const TimeNodePDE &tn20) const;
    void border2_calculate(DoubleMatrix &u20, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn20) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    void gridMethod0(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod1(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
    void gridMethod2(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep);
};

#endif // HYPERBOLICIBVP_H
