#ifndef HYPERBOLIC_IBVP_H
#define HYPERBOLIC_IBVP_H

#include "ibvp.h"

class MINIMUMSHARED_EXPORT IHyperbolicIBVP : public InitialBoundaryValueProblemPDE
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IFinalIHyperbolicIBVP : public FinalBoundaryValueProblemPDE
{
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IWaveEquationIBVP : public IHyperbolicIBVP
{
public:
    explicit IWaveEquationIBVP(double waveSpeed = 1.0, double waveDissipation = 0.0);
    virtual ~IWaveEquationIBVP();
    IWaveEquationIBVP(const IWaveEquationIBVP &);
    IWaveEquationIBVP & operator =(const IWaveEquationIBVP &);

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    virtual void setWaveSpeed(double waveSpeed);
    virtual void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

protected:
    virtual double lambda() const;
    double _waveSpeed;
    double _waveDissipation;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IFinalWaveEquationIBVP : public IFinalIHyperbolicIBVP
{
public:
    explicit IFinalWaveEquationIBVP(double waveSpeed = 1.0, double waveDissipation = 0.0);
    IFinalWaveEquationIBVP(const IFinalWaveEquationIBVP &);
    IFinalWaveEquationIBVP & operator =(const IFinalWaveEquationIBVP &);
    virtual ~IFinalWaveEquationIBVP();

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    virtual void setWaveSpeed(double waveSpeed);
    virtual void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

    virtual double lambda() const;

protected:
    double _waveSpeed;
    double _waveDissipation;
};

//--------------------------------------------------------------------------------------------------------------//

#endif // HYPERBOLIC_IBVP_H
