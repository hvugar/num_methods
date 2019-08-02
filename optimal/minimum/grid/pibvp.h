#ifndef PARABOLIC_IBVP_H
#define PARABOLIC_IBVP_H

#include "ibvp.h"

class MINIMUMSHARED_EXPORT IParabolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IFinalParabolicIBVP : public FinalBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IHeatEquationIBVP : public IParabolicIBVP
{
public:
    explicit IHeatEquationIBVP(double thermalDiffusivity = 1.0);
    virtual ~IHeatEquationIBVP();
    IHeatEquationIBVP(const IHeatEquationIBVP &);
    IHeatEquationIBVP & operator =(const IHeatEquationIBVP &);

    virtual double thermalDiffusivity() const;
    virtual void setThermalDiffusivity(double thermalDiffusivity);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

    virtual double env0(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return 0.0; }
    virtual double env1(const SpaceNodePDE &sn, const TimeNodePDE &tn) const { return 0.0; }

    void calculateU(DoubleMatrix &u, double a, double alpha, double lambda);
    void gridMethod(DoubleVector &u) const;
    void gridMethod(DoubleVector &u, double a = 1.0) const;
    void calculateMVD(DoubleMatrix &u) const;
    /* dirichlet conditions */
    void calculateN2L2RD(DoubleMatrix &u) const;
    void calculateN4L2RD(DoubleMatrix &u) const;
    void calculateN4L2RDX(DoubleMatrix &u) const;
    void calculateN6L2RD(DoubleMatrix &u) const;

protected:
    virtual double lambda() const;
    double _thermalDiffusivity; // температуропроводность
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IFinalHeatEquationIBVP : public IFinalParabolicIBVP
{
protected:
    explicit IFinalHeatEquationIBVP(double thermalDiffusivity = 1.0);
    IFinalHeatEquationIBVP(const IFinalHeatEquationIBVP &);
    IFinalHeatEquationIBVP & operator =(const IFinalHeatEquationIBVP &);
    virtual ~IFinalHeatEquationIBVP();

    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double thermalDiffusivity() const;
    virtual void setThermalDiffusivity(double thermalDiffusivity);

    void gridMethod(DoubleVector &u, SweepMethodDirection direction = ForwardSweep) const;
    void gridMethod(DoubleMatrix &u, SweepMethodDirection direction = ForwardSweep) const;

    void calculate(DoubleCube &u, SweepMethodDirection direction = ForwardSweep) const;

protected:
    virtual double lambda() const;
    double _thermalDiffusivity; // Температуропроводность
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

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

#endif // PARABOLIC_IBVP_H

// a - thermal diffusivity
// k - thermal conductivity
// c - specific heat capacity
// p - density
