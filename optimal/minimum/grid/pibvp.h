#ifndef PARABOLIC_IBVP_H
#define PARABOLIC_IBVP_H

#include "ibvp.h"
#include "../deltagrid.h"
#include <algorithm>
#include <thread>
#include <iostream>

#define PARABOLIC_IBVP_H_D2V1_BORDER_O2
#define PARABOLIC_IBVP_H_D2V1_FX_X
#define PARABOLIC_IBVP_H_D2V1_FX_Y
//#define PARABOLIC_IBVP_H_D2V1_BR_X
#define PARABOLIC_IBVP_H_D2V1_BR_Y

#define PARABOLIC_IBVP_H_D1V1_BORDER_O2
//#define PARABOLIC_IBVP_H_D1V1_FX
//#define PARABOLIC_IBVP_H_D1V1_BR_LEFT
//#define PARABOLIC_IBVP_H_D1V1_BR_RIGHT


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

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const = 0;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const = 0;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const = 0;

    virtual double weight() const;  //

    double _thermalDiffusivity;     // температуропроводность
    double _thermalConductivity;    // теплопроводность - heat transfer coefficient
    double _thermalConvection;      //
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
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

struct MINIMUMSHARED_EXPORT LoadedSpacePoint : public SpacePoint
{
    LoadedSpacePoint(double x = 0.0, double y = 0.0, double z = 0.0, double d = 0.0) : SpacePoint(x, y, z), d(d) {}
    double d;
};

class MINIMUMSHARED_EXPORT ILoadedHeatEquationIBVP : public IHeatEquationIBVP
{
public:
    explicit ILoadedHeatEquationIBVP(double thermalDiffusivity = 1.0, double thermalConductivity = 0.0, double thermalConvection = 0.0);
    virtual ~ILoadedHeatEquationIBVP();
    ILoadedHeatEquationIBVP(const ILoadedHeatEquationIBVP &);
    ILoadedHeatEquationIBVP& operator=(const ILoadedHeatEquationIBVP &);

    virtual void explicit_calculate_D1V1() const;// TO-DO
    virtual void implicit_calculate_D1V1() const;// TO-DO

    virtual void explicit_calculate_D2V1() const;// TO-DO
    virtual void implicit_calculate_D2V1() const;// TO-DO

    void setLoadedPoints(const std::vector<LoadedSpacePoint> &loadedPoints);
    const std::vector<LoadedSpacePoint> loadedPoints() const;

private:
    std::vector<LoadedSpacePoint> _loadedPoints;
};

//------------------------------------------------------------------------------------------------------------------------------------------------------//

class MINIMUMSHARED_EXPORT ILoadedHeatEquationFBVP : public IHeatEquationFBVP
{
public:
    explicit ILoadedHeatEquationFBVP(double thermalDiffusivity = 1.0, double thermalConductivity = 0.0, double thermalConvection = 0.0);
    ILoadedHeatEquationFBVP(const ILoadedHeatEquationFBVP &);
    ILoadedHeatEquationFBVP& operator=(const ILoadedHeatEquationFBVP &);
    virtual ~ILoadedHeatEquationFBVP();

    virtual void explicit_calculate_D1V1() const;// TO-DO
    virtual void implicit_calculate_D1V1() const;// TO-DO

    virtual void explicit_calculate_D2V1() const;// TO-DO
    virtual void implicit_calculate_D2V1() const;// TO-DO

    void implicit_calculate_D2V2() const;// TO-DO

    void setLoadedPoints(const std::vector<LoadedSpacePoint> &loadedPoints);
    const std::vector<LoadedSpacePoint> loadedPoints() const;

private:
    std::vector<LoadedSpacePoint> _loadedPoints;

};


#endif // PARABOLIC_IBVP_H

// a - thermal diffusivity
// k - thermal conductivity
// c - specific heat capacity
// p - density
