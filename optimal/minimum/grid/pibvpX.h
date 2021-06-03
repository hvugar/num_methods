#ifndef PARABOLIC_IBVPX_H
#define PARABOLIC_IBVPX_H

#include "pibvp.h"

enum SweepMethodDirection
{
    ForwardSweep = 1,
    BackwardSweep = 2,
    Centered = 3
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IHeatEquationIBVPEx : public IHeatEquationIBVP
{
public:
    virtual double env0(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }
    virtual double env1(const SpaceNodePDE &, const TimeNodePDE &) const { return 0.0; }

    void calculateU(DoubleMatrix &u, double a, double alpha, double weight);
    //void gridMethod(DoubleVector &u) const;
    void gridMethod(DoubleVector &u, double a = 1.0) const;
    void calculateMVD(DoubleMatrix &u) const;
    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u) const;
    void calculateN4L2RD(DoubleMatrix &u) const;
    void calculateN4L2RDX(DoubleMatrix &u) const;
    void calculateN6L2RD(DoubleMatrix &u) const;
};

class MINIMUMSHARED_EXPORT IHeatEquationFBVPEx : public IHeatEquationFBVP
{
public:
    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;
    void calculate(DoubleCube &u, SweepMethodDirection direction = ForwardSweep) const;
};

class MINIMUMSHARED_EXPORT NewtonHeatEquation : public IParabolicIBVP
{
public:
    double lambda0;
    double lambda1;
    double lambda2;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
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
};

#endif // PARABOLIC_IBVPX_H
