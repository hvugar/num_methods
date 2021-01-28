#ifndef PARABOLIC_IBVP_H
#define PARABOLIC_IBVP_H

#include "ibvp.h"
#include "../deltagrid.h"
#define PARABOLIC_IBVP_O2

/**
 * @brief The IParabolicIBVP class
 * @class IParabolicIBVP
 * @see InitialBoundaryValueProblemPDE
 */
class MINIMUMSHARED_EXPORT IParabolicIBVP : public InitialBoundaryValueProblemPDE
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
 * @brief The IParabolicFBVP class
 * @class IParabolicFBVP
 * @see FinalBoundaryValueProblemPDE
 */
class MINIMUMSHARED_EXPORT IParabolicFBVP : public FinalBoundaryValueProblemPDE
{
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual void layerInfo(const DoubleVector &, const TimeNodePDE &) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IHeatEquationIBVP : public IParabolicIBVP
{
public:
    explicit IHeatEquationIBVP(double thermalDiffusivity = 1.0, double thermalConductivity = 0.0, double thermalConvection = 0.0);
    virtual ~IHeatEquationIBVP();
    IHeatEquationIBVP(const IHeatEquationIBVP &);
    IHeatEquationIBVP& operator=(const IHeatEquationIBVP &);

    virtual double thermalDiffusivity() const;
    virtual void setThermalDiffusivity(double thermalDiffusivity);

    virtual double thermalConvection() const;
    virtual void setThermalConvection(double thermalConvection);

    virtual double thermalConductivity() const;
    virtual void setThermalConductivity(double thermalConductivity);

    virtual void explicit_calculate_D1V1() const;// TO-DO
    virtual void implicit_calculate_D1V1() const;// COMPLETED

    virtual void explicit_calculate_D2V1() const;// TO-DO
    virtual void implicit_calculate_D2V1() const;// COMPLETED

    void init() const;
    void next1() const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double weight() const;  //
    double _thermalDiffusivity;     // температуропроводность
    double _thermalConductivity;    // теплопроводность - heat transfer coefficient
    double _thermalConvection;      //

public:
    bool _userHalfValues = false;

private:
    void *initParams;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT IHeatEquationFBVP : public IParabolicFBVP
{
public:
    explicit IHeatEquationFBVP(double thermalDiffusivity = 1.0, double thermalConductivity = 0.0, double thermalConvection = 0.0);
    IHeatEquationFBVP(const IHeatEquationFBVP &);
    IHeatEquationFBVP& operator=(const IHeatEquationFBVP &);
    virtual ~IHeatEquationFBVP();

    virtual double thermalDiffusivity() const;
    virtual void setThermalDiffusivity(double thermalDiffusivity);

    virtual double thermalConvection() const;
    virtual void setThermalConvection(double thermalConvection);

    virtual double thermalConductivity() const;
    virtual void setThermalConductivity(double thermalConductivity);

    virtual void explicit_calculate_D1V1() const;// TO-DO
    virtual void implicit_calculate_D1V1() const;// COMPLETED

    virtual void explicit_calculate_D2V1() const;// TO-DO
    virtual void implicit_calculate_D2V1() const;// COMPLETED

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double weight() const;  //
    double _thermalDiffusivity;     // температуропроводность
    double _thermalConductivity;    // теплопроводность - heat transfer coefficient
    double _thermalConvection;      //

public:
    bool _userHalfValues = false;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT ILoadedHeatEquationIBVP : public IHeatEquationIBVP
{
public:
    explicit ILoadedHeatEquationIBVP(double thermalDiffusivity = 1.0, double thermalConductivity = 0.0, double thermalConvection = 0.0);
    virtual ~ILoadedHeatEquationIBVP();
    ILoadedHeatEquationIBVP(const ILoadedHeatEquationIBVP &);
    ILoadedHeatEquationIBVP& operator=(const ILoadedHeatEquationIBVP &);

    virtual void explicit_calculate_D1V1() const;// TO-DO
    virtual void implicit_calculate_D1V1() const;// COMPLETED

    virtual void explicit_calculate_D2V1() const;// TO-DO
    virtual void implicit_calculate_D2V1() const;// COMPLETED

    void setLoadedPoints(const std::vector<SpacePoint> &loadedPoints);
    const std::vector<SpacePoint> loadedPoints() const;
    std::vector<double> _d;

private:
    std::vector<SpacePoint> _loadedPoints;
};


#endif // PARABOLIC_IBVP_H

// a - thermal diffusivity
// k - thermal conductivity
// c - specific heat capacity
// p - density
