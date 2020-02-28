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
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

/**
 * @brief The IFinalIHyperbolicIBVP class
 * @class IFinalIHyperbolicIBVP
 * @see FinalBoundaryValueProblemPDE
 */
class MINIMUMSHARED_EXPORT IHyperbolicFBVP : public FinalBoundaryValueProblemPDE
{
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IWaveEquationIBVP : public IHyperbolicIBVP
{
public:
    explicit IWaveEquationIBVP(double waveSpeed = 1.0, double waveDissipation = 0.0, double unknownB = 0.0, double restoration = 0.0);
    virtual ~IWaveEquationIBVP();
    IWaveEquationIBVP(const IWaveEquationIBVP &);
    IWaveEquationIBVP& operator=(const IWaveEquationIBVP &);

    double waveSpeed() const;
    double waveDissipation() const;

    void setWaveSpeed(double waveSpeed);
    void setWaveDissipation(double waveDissipation);

    void setUnknownB(double unknownB);
    double unknownB() const;

    void setRestoration(double restoration);
    double restoration() const;

    void explicit_calculate_D1V1() const; // TO-DO
    void implicit_calculate_D1V1() const; // TO-DO

    void explicit_calculate_D2V1() const; // TO-DO
    void implicit_calculate_D2V1() const; // TO-DO

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double weight() const;
    double _waveSpeed;
    double _waveDissipation;
    double _unknownB;
    double _restoration;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IWaveEquationFBVP : public IHyperbolicFBVP
{
public:
    explicit IWaveEquationFBVP(double waveSpeed = 1.0, double waveDissipation = 0.0, double unknownB = 0.0, double restoration = 0.0);
    IWaveEquationFBVP(const IWaveEquationFBVP &);
    IWaveEquationFBVP& operator=(const IWaveEquationFBVP &);
    virtual ~IWaveEquationFBVP();

    virtual double waveSpeed() const;
    virtual double waveDissipation() const;

    virtual void setWaveSpeed(double waveSpeed);
    virtual void setWaveDissipation(double waveDissipation);

    void setUnknownB(double unknownB);
    double unknownB() const;

    void setRestoration(double restoration);
    double restoration() const;

    void explicit_calculate_D1V1() const; // TO-DO
    void implicit_calculate_D1V1() const; // TO-DO

    void explicit_calculate_D2V1() const; // TO-DO
    void implicit_calculate_D2V1() const; // TO-DO

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double weight() const;
    double _waveSpeed;
    double _waveDissipation;
    double _unknownB;
    double _restoration;
};

//--------------------------------------------------------------------------------------------------------------//

#endif // HYPERBOLIC_IBVP_H
