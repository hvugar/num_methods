#ifndef SOLVER5_H
#define SOLVER5_H

#include "global.h"

namespace p3p5
{

class Functional;

class PROBLEM3P_SHARED_EXPORT LoadedHeatEquationIBVP : public ILoadedHeatEquationIBVP
{
public:
    LoadedHeatEquationIBVP(Functional*);
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

private:
    double _initial_temperature = 0.0;
    double _enviroment_temperature = 0.5;

    Functional* _functional;
    double _fx;
};

class PROBLEM3P_SHARED_EXPORT Functional /*: public RnFunction, public IGradient, public IProjection, public IPrinter*/
{
public:
    static void Main(int argc, char** argv);

    double lambda0 = 0.10;
    double lambda1 = 0.01;
    double initial_temperature = 0.0;
    double envrmnt_temperature = 1.0;

    size_t heating_point_number = 2;

    auto q(const TimeNodePDE &tn, size_t i) const -> double;
    auto z(const TimeNodePDE &tn, size_t i) const -> SpacePoint;

    auto timeDimension() const -> Dimension { return Dimension(0.01, 0, 1000); }

    auto spaceDimensionX() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionY() const -> Dimension { return Dimension(0.01, 0, 100); }

    auto spaceDimensionZ() const -> Dimension { return Dimension(0.01, 0, 100); }
};



}



#endif // SOLVER5_H
