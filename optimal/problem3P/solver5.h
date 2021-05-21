#ifndef SOLVER5_H
#define SOLVER5_H

#include "global.h"

namespace p3p5
{

class Functional;

class Shared
{
public:
    double lambda1 = 0.01;
    double initial_temperature = 0.0;
    double envrmnt_temperature = 5.0;

    size_t heating_point_number = 2;

    DoubleMatrix U;
    DoubleMatrix V = DoubleMatrix(100, 100, 0.0);

    auto q(const TimeNodePDE &tn, size_t i) const -> double;
    auto z(const TimeNodePDE &tn, size_t i) const -> SpacePoint;

protected:

    Dimension _timeDimension   = Dimension(0.01, 0, 1000);
    Dimension _spaceDimensionX = Dimension(0.01, 0, 100);
    Dimension _spaceDimensionY = Dimension(0.01, 0, 100);
    Dimension _spaceDimensionZ = Dimension(0.01, 0, 100);
};

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationIBVP : public IHeatEquationIBVP, public virtual Shared /*public ILoadedHeatEquationIBVP*/
{
public:
    LoadedHeatEquationIBVP(double a = 1.0, double lambda0 = 0.0, double lambda1 = 0.0);
    LoadedHeatEquationIBVP(const LoadedHeatEquationIBVP &);
    LoadedHeatEquationIBVP & operator =(const LoadedHeatEquationIBVP &);
    virtual ~LoadedHeatEquationIBVP() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition ic) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;
};

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationFBVP : public IHeatEquationFBVP, public virtual Shared /*public ILoadedHeatEquationFBVP*/
{
public:
    LoadedHeatEquationFBVP(double a = -1.0, double lambda0 = 0.0, double lambda1 = 0.0);
    LoadedHeatEquationFBVP(const LoadedHeatEquationFBVP &);
    LoadedHeatEquationFBVP & operator =(const LoadedHeatEquationFBVP &);
    virtual ~LoadedHeatEquationFBVP() override;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition fc) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &bc) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

    virtual auto layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;
};


class PROBLEM3P_SHARED_EXPORT Functional : public LoadedHeatEquationIBVP, public LoadedHeatEquationFBVP /*: public RnFunction, public IGradient, public IProjection, public IPrinter*/
{
public:
    static void Main(int argc, char** argv);

    Functional(double diffusivity, double convection, double conductivity = 0.0, double lambda = 0.0);
};

}



#endif // SOLVER5_H
