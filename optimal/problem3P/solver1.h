#ifndef PROBLEM3P_SOLVER1_H
#define PROBLEM3P_SOLVER1_H

#include "global.h"

namespace p3p1
{

struct HeatSourceParams
{
    double *v1 = nullptr;
    double *v2 = nullptr;
    double *q1 = nullptr;
    double *q2 = nullptr;
};

class Solver1;

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP, virtual public ISecondOrderLinearODEIBVP
{
public:
    HeatEquationIBVP(Solver1 *solver = nullptr);
    virtual ~HeatEquationIBVP() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;
    virtual auto weight() const -> double override;

public:
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    virtual auto A(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto B(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto C(const PointNodeODE &node, size_t row) const -> double override;
    virtual auto D(const PointNodeODE &node, size_t row) const -> double;

protected:
    virtual auto initial(InitialCondition condition, size_t row = 1) const  -> double override;
    virtual auto boundary(const PointNodeODE &, BoundaryConditionPDE &, size_t) const -> double override { return 0.0; }
    virtual auto count() const  -> size_t override;

public:
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

    size_t i;

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

    void setPointNumber(size_t heatSourceNumber, size_t measrPointNumber);

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    //virtual void frw_calculate() const;

    size_t heatSourceNumber = 2;
    size_t measrPointNumber = 4;
    SpacePoint *measurePoints = new SpacePoint[measrPointNumber];
    SpacePoint *measurePointValues = new SpacePoint[measrPointNumber];
    double lambda0 = 0.001;

    DoubleMatrix alpha1;
    DoubleMatrix alpha2;
    DoubleMatrix alpha3;

    DoubleMatrix betta1;
    DoubleMatrix betta2;
    DoubleMatrix betta3;

    DoubleMatrix nominU;

    double **q = nullptr;
    double **v = nullptr;
    SpacePoint **z = nullptr;

    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;

private:
    double _initialValue = 0.0;
    double _environmentTemperature = 0.0;
};

}

#endif // PROBLEM3P_SOLVER1_H
