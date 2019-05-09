#ifndef PARABOLIC_IBVP_H
#define PARABOLIC_IBVP_H

#include "ibvp.h"
#include "../printer.h"

// a - thermal diffusivity
// k - thermal conductivity
// c - specific heat capacity
// p - density

class MINIMUMSHARED_EXPORT IParabolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
};

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IHeatEquationIBVP : public IParabolicIBVP
{
public:
    IHeatEquationIBVP(double thermalDiffusivity = 1.0);

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) {}

    void gridMethod(DoubleVector &u, double a = 1.0, SweepMethodDirection direction = ForwardSweep);
    void gridMethod1(DoubleVector &u, double a = 1.0);

protected:
    double _thermalDiffusivity;
};

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT CCIParabolicIBVP : public IParabolicIBVP
{
protected:
    virtual void layerInfo(const DoubleVector &, unsigned int) const;
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const;

public:
    void calculate1(DoubleVector &u, double a = 1.0) const;
    void calculate1(DoubleMatrix &u, double a = 1.0) const;
};

/**
 * @brief The ParabolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t), t in (0,T], x in (0, l),
 * u(x,0) = fi(x), x in [0,l],
 * u(0,t) = m1(t), t in (0,T],
 * u(l,t) = m2(t), t in (0,T].
 */
class MINIMUMSHARED_EXPORT ParabolicIBVP : public IParabolicIBVP
{
protected:
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:

    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;

    void gridMethod1L(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod1LT(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;

    void gridMethod1R(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod11(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod2(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;

    void calculateMVD(DoubleMatrix &u) const;
    void calculateMVD_TEST(DoubleMatrix &u) const;

    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u) const;
    //void calculateN2R2LD(DoubleMatrix &u) const;

    void calculateN4L2RD(DoubleMatrix &u) const;
    void calculateN4L2RDX(DoubleMatrix &u) const;

    //void calculateN4R2LD(DoubleMatrix &u) const;

    void calculateN6L2RD(DoubleMatrix &u) const;
    //void calculateN6R2LD(DoubleMatrix &u) const;
};

#endif // PARABOLIC_IBVP_H
