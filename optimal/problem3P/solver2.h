#ifndef SOLVER2_H
#define SOLVER2_H


#include "global.h"

namespace p3p0
{

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP
{
public:
    HeatEquationIBVP();
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

public:
    auto q(size_t i, const TimeNodePDE &tn) const -> double;
    auto v(size_t i, const PointNodeODE &tn, SpacePoint &vl) const -> void;

public:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

public:
    virtual auto layerInfo(const DoubleVector &, const TimeNodePDE &) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

public:
    size_t i;

    double lambda1 = +0.01;
    double theta = +2.0;
};

class PROBLEM3P_SHARED_EXPORT Solver2
{
public:
    static void Main(int argc, char** argv);

public:
    Solver2();
};

}

#endif // SOLVER2_H
