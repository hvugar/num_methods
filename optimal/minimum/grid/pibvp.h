#ifndef PARABOLIC_IBVP_H
#define PARABOLIC_IBVP_H

#include "ibvp.h"

// a - thermal diffusivity
// k - thermal conductivity
// c - specific heat capacity
// p - density

class MINIMUMSHARED_EXPORT IParabolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
};

//--------------------------------------------------------------------------------------------------------------//

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
    void gridMethod(DoubleVector &u) const;
    void gridMethod(DoubleVector &u, double a = 1.0) const;
    void calculateMVD(DoubleMatrix &u) const;

    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u) const;
    void calculateN4L2RD(DoubleMatrix &u) const;
    void calculateN6L2RD(DoubleMatrix &u) const;
    void calculateN4L2RDX(DoubleMatrix &u) const;
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

    virtual double env0(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double env1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

public:
    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}



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


protected:
    double _thermalDiffusivity;
};

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT NewtonHeatEquation : public IParabolicIBVP
{
public:
    double lambda0;
    double lambda1;
    double lambda2;

protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
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

//--------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IBackwardParabolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
};

/**
 * @brief The BackwardParabolicIBVP class
 * u_t(x,t) = a(x,t)u_xx(x,t) + f(x,t),
 *          t in [0,T), x in (0, l),
 * u(x,T) = fi(x), x in [0,l],
 * u(0,t) = m1(t), t in [0,T),
 * u(l,t) = m2(t), t in [0,T).
 */
class MINIMUMSHARED_EXPORT BackwardParabolicIBVP : public IBackwardParabolicIBVP
{
protected:
    virtual double initial(const SpaceNodePDE &sn) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;
    virtual double a(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, unsigned int) const {}
    virtual void layerInfo(const DoubleMatrix &, unsigned int) const {}

public:

    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;

    void calculate(DoubleCube &u, SweepMethodDirection direction = ForwardSweep) const;
};

#endif // PARABOLIC_IBVP_H
