#ifndef PROBLEM3P_SOLVER1_H
#define PROBLEM3P_SOLVER1_H

#include "global.h"

namespace p3p
{

class Solver1;

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP, virtual public IFirstOrderLinearODEIVP
{
public:
    HeatEquationIBVP(Solver1 *solver = nullptr);
    virtual ~HeatEquationIBVP();

    virtual auto initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto weight() const -> double override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const  -> double override;
    virtual auto B(const PointNodeODE &node, unsigned int row = 1) const  -> double override;
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double;
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const  -> double override;
    virtual auto count() const  -> unsigned int override;
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

    unsigned int i;
    DoubleMatrix lastLayerU;

public:
    Solver1 *solver;
};

class Solver1
{
public:
    static void Main(int argc, char** argv);

public:
    Solver1();
    virtual ~Solver1();

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void frw_calculate() const;

    const size_t heatSourceNumber = 2;
    std::vector<SpacePoint*> heatSourceRoutes;
    const size_t measurePointNumber = 4;
    SpacePoint *measurePoints = new SpacePoint[measurePointNumber];
    SpacePoint *measurePointValues = new SpacePoint[measurePointNumber];
    double environmentTemperature = 0.0;
    double lambda0 = 0.001;

    double *alpha1 = new double[heatSourceNumber*measurePointNumber];
    double *alpha2 = new double[heatSourceNumber*measurePointNumber];
    double *alpha3 = new double[heatSourceNumber*measurePointNumber];
    double *betta1 = new double[heatSourceNumber*measurePointNumber];
    double *betta2 = new double[heatSourceNumber*measurePointNumber];
    double *betta3 = new double[heatSourceNumber*measurePointNumber];
};

}

#endif // PROBLEM3P_SOLVER1_H
