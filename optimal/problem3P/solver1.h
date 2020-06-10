#ifndef PROBLEM3P_SOLVER1_H
#define PROBLEM3P_SOLVER1_H

#include "global.h"

namespace p3p1
{

struct ProblemParams
{
    double *q;
    double *v;
    SpacePoint *z;
    SpacePoint *zt;
    SpacePoint *u;
    SpacePoint *p;
    SpacePoint *f;

};

class Solver1;

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP, virtual public ISecondOrderLinearODEIVP
{
public:
    HeatEquationIBVP(Solver1 *solver = nullptr);
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

protected:
    virtual auto initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

public:
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    virtual auto A(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto B(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto C(const PointNodeODE &node, size_t row) const -> double override;

protected:
    virtual auto initial(InitialCondition condition, size_t row = 1) const  -> double override;
    virtual auto count() const  -> size_t override;

public:
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

public:
    Solver1 *solver;
    size_t i;
};

class PROBLEM3P_SHARED_EXPORT HeatEquationFBVP : virtual public IHeatEquationFBVP, virtual public ISecondOrderLinearODEFVP
{
public:
    HeatEquationFBVP(Solver1 *solver = nullptr);
    HeatEquationFBVP(const HeatEquationFBVP &);
    HeatEquationFBVP & operator =(const HeatEquationFBVP &);
    virtual ~HeatEquationFBVP() override;

protected:
    virtual auto final(const SpaceNodePDE &sn, FinalCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;

public:
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    virtual auto A(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto B(const PointNodeODE &node, size_t row, size_t col) const -> double override;
    virtual auto C(const PointNodeODE &node, size_t row) const -> double override;

protected:
    virtual auto final(FinalCondition condition, size_t row = 1) const  -> double override;
    virtual auto count() const  -> size_t override;

public:
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

public:
    Solver1 *solver;
    size_t i;
};

class PROBLEM3P_SHARED_EXPORT Solver1
{
public:
    static void Main(int argc, char** argv);

public:
    Solver1(const Dimension &timeDimension,
            const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);
    virtual ~Solver1();

    void setPointNumber(size_t heatSourceNumber, size_t measrPointNumber);

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;

    virtual double A1(const PointNodeODE &node, size_t row, size_t col, size_t i) const;
    virtual double A2(const PointNodeODE &node, size_t row, size_t col, size_t i) const;
    virtual double A3(const PointNodeODE &node, size_t row, size_t i) const;
    virtual double A4(const PointNodeODE &node, size_t row, size_t i) const;

    size_t heatSourceNumber = 2;
    size_t measrPointNumber = 4;
    double lambda0 = 0.001;
    double lambda = 0.01;

    DoubleMatrix alpha1;
    DoubleMatrix alpha2;
    DoubleMatrix alpha3;
    DoubleMatrix betta1;
    DoubleMatrix betta2;
    DoubleMatrix betta3;
    DoubleMatrix nU;
    SpacePoint *measurePoints;

    ProblemParams *sourceParams = nullptr;

    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }

    DoubleMatrix V, U;

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;

private:
    double frw_initialValue = 1.0;
    double environmentTemperature = 1.0;
};

}

#endif // PROBLEM3P_SOLVER1_H
