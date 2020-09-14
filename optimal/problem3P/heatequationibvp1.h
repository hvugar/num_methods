#ifndef HEATEQUATIONIBVP1_H
#define HEATEQUATIONIBVP1_H

#include "global.h"

class HeatEquationIBVP1 : public IHeatEquationIBVP
{
public:
    HeatEquationIBVP1();
    HeatEquationIBVP1(const HeatEquationIBVP1 &);
    HeatEquationIBVP1 & operator =(const HeatEquationIBVP1 &);
    virtual ~HeatEquationIBVP1() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;

protected:
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

protected:
    //virtual auto A(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    //virtual auto B(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    //virtual auto C(const PointNodeODE &node, size_t row) const -> double override;
    //virtual auto initial(InitialCondition condition, size_t row = 1) const  -> double override;
    //virtual auto count() const  -> size_t override;
    //virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;
    //virtual auto dimension() const -> Dimension override;

public:
    auto q(size_t i, const TimeNodePDE &tn) const -> double;
    auto v(size_t i, const PointNodeODE &tn) const -> double;
    auto theta(const TimeNodePDE &tn) const -> double;
    auto frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const -> void;

public:
    size_t i;
    double lambda0 = 0.0;
    double lambda1 = +0.01;
};

#endif // HEATEQUATIONIBVP1_H
