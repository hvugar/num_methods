#ifndef HYPERBOLIC_IBVP_H
#define HYPERBOLIC_IBVP_H

#include "ibvp.h"

/**
 * @brief The IHyperbolicIBVP class
 * @class IHyperbolicIBVP
 * @see InitialBoundaryValueProblemPDE
 */
class MINIMUMSHARED_EXPORT IHyperbolicIBVP : public InitialBoundaryValueProblemPDE
{
public:
    IHyperbolicIBVP();
    IHyperbolicIBVP(const IHyperbolicIBVP &);
    IHyperbolicIBVP & operator = (const IHyperbolicIBVP &);
    virtual ~IHyperbolicIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const {}
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const {}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

/**
 * @brief The IFinalIHyperbolicIBVP class
 * @class
 * @see
 */
class MINIMUMSHARED_EXPORT IHyperbolicFBVP : public FinalBoundaryValueProblemPDE
{
public:
    IHyperbolicFBVP();
    IHyperbolicFBVP(const IHyperbolicFBVP &);
    IHyperbolicFBVP & operator = (const IHyperbolicFBVP &);
    virtual ~IHyperbolicFBVP();

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

    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    virtual void setWaveSpeed(double waveSpeed);
    virtual void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

protected:
    virtual double weight() const;
    double _waveSpeed;
    double _waveDissipation;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IWaveEquationFBVP : public IHyperbolicFBVP
{
public:
    explicit IWaveEquationFBVP(double waveSpeed = 1.0, double waveDissipation = 0.0);
    IWaveEquationFBVP(const IWaveEquationFBVP &);
    IWaveEquationFBVP & operator =(const IWaveEquationFBVP &);
    virtual ~IWaveEquationFBVP();

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    virtual void setWaveSpeed(double waveSpeed);
    virtual void setWaveDissipation(double waveDissipation);

    void explicit_calculate_D1V1() const;
    void implicit_calculate_D1V1() const;

    void explicit_calculate_D2V1() const;
    void implicit_calculate_D2V1() const;

protected:
    virtual double weight() const;
    double _waveSpeed;
    double _waveDissipation;
};

//--------------------------------------------------------------------------------------------------------------//

#endif // HYPERBOLIC_IBVP_H
