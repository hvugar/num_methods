#ifndef PROBLEM3P_SOLVER1_H
#define PROBLEM3P_SOLVER1_H

#include "global.h"
#include <functional>   // std::bind

namespace p3p1
{

struct SpacePointX : SpacePoint
{
    SpacePointX(double x = .0, double y = .0, double dx = .0, double dy = .0)
        : SpacePoint(x, y), dx(dx), dy(dy) {}
    SpacePointX(DoubleVector &z) : SpacePoint(z[0], z[1]), dx(z[2]), dy(z[3]) {}

    void toDoubleVector(DoubleVector &z) { z[0] = x; z[1] = y; z[2] = dx; z[3] = dy; }

    double dx = .0;
    double dy = .0;
};

struct ProblemParams
{
    double *q;  // istilik menbeyinin gucu
    double *vX; // suert
    double *vY; // suert

    SpacePointX *z;
    SpacePointX *f;

    SpacePoint *u;
    SpacePoint *p;
};

class Solver1;

//class HeatEquationIBVP2 : virtual public IHeatEquationIBVP, virtual public ISecondOrderLinearODEIVP
//{
//public:
//    HeatEquationIBVP2(Solver1 *solver = nullptr);
//    HeatEquationIBVP2(const HeatEquationIBVP2 &);
//    HeatEquationIBVP2 & operator =(const HeatEquationIBVP2 &);
//    virtual ~HeatEquationIBVP2() override;

//public:
//    virtual auto initial(const SpaceNodePDE &/*sn*/, InitialCondition /*condition*/) const -> double override { return 0.5; }
//    virtual auto boundary(const SpaceNodePDE &/*sn*/, const TimeNodePDE &/*tn*/, BoundaryConditionPDE &cn) const -> double override
//    { cn = BoundaryConditionPDE::Robin(_lambda, -1.0, _lambda); return 0.5; }
//    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override
//    {
//        size_t ln = static_cast<size_t>(tn.i);

//        const ProblemParams &pp = sourceParams[ln];
//        double fx = _lambda0 * _environmentTemperature;

//        double sum = 0.0;
//        for (size_t i=0; i<heatSourceNumber; i++)
//        {
//            const SpacePointX &zi = pp.z[i];
//            if ( isPointOnPlate(zi) )
//            {
//                sum += pp.q[i] * DeltaFunction::gaussian(sn, zi, SpacePoint(spaceDimensionX().step()*_factor, spaceDimensionY().step()*_factor));
//            }
//            else
//            {
//                std::cout << "frw_f" << tn.t << ", " << zi.x << ", " << zi.y << std::endl;
//            }
//        }

//        return fx + sum;
//    }

//public:
//    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
//    virtual auto timeDimension() const -> Dimension override;
//    virtual auto spaceDimensionX() const -> Dimension override;
//    virtual auto spaceDimensionY() const -> Dimension override;
//    virtual auto spaceDimensionZ() const -> Dimension override;

//    virtual auto A(const PointNodeODE &node, size_t row, size_t col) const -> double override;
//    virtual auto B(const PointNodeODE &node, size_t row, size_t col) const -> double override;
//    virtual auto C(const PointNodeODE &node, size_t row) const -> double override;

//protected:
//    virtual auto initial(InitialCondition condition, size_t row = 1) const  -> double override;
//    virtual auto count() const  -> size_t override;

//public:
//    virtual auto dimension() const -> Dimension override;
//    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

//public:
//    size_t heatSourceNumber = 2;
//    size_t measrPointNumber = 4;
//    double _lambda0 = +0.001;
//    double _lambda1 = -0.01;
//    double _factor = 1.0;

//    size_t i;
//};

class HeatEquationIBVP : virtual public IHeatEquationIBVP, virtual public ISecondOrderLinearODEIVP
{
public:
    HeatEquationIBVP(Solver1 *solver = nullptr);
    HeatEquationIBVP(const HeatEquationIBVP &);
    HeatEquationIBVP & operator =(const HeatEquationIBVP &);
    virtual ~HeatEquationIBVP() override;

public:
    auto q(size_t i, const TimeNodePDE &tn) const -> double;
    auto v(size_t i, const PointNodeODE &tn) const -> double;

public:
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
    bool findUMode = true;
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
    //virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

public:
    Solver1 *solver;
    size_t i;
};

class PROBLEM3P_SHARED_EXPORT Solver1 : public IGradient, public RnFunction, public IProjection, public IPrinter
{
public:
    static void Main(int argc, char** argv);
    static void optimize(int argc, char **argv);

public:
    Solver1(const Dimension &timeDimension,
            const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);
    virtual ~Solver1();

    void setPointNumber(size_t heatSourceNumber, size_t measrPointNumber);

    virtual void gradient(const DoubleVector &x, DoubleVector &g) const;
    virtual double fx(const DoubleVector &x) const;
    auto integral(const DoubleMatrix &) const -> double;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &x) const -> void;
    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g,
                       double f, double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const -> double;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void bcw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    virtual double A1(const PointNodeODE &node, size_t row, size_t col, size_t i) const;
    virtual double A2(const PointNodeODE &node, size_t row, size_t col, size_t i) const;
    virtual double A3(const PointNodeODE &node, size_t row, size_t i) const;
    virtual double A4(const PointNodeODE &node, size_t row, size_t col, size_t i) const;

    bool isPointOnPlate(const SpacePoint &z) const;

    void drawImages(Solver1 &solver) const;
    double minU, maxU;
    bool saveMinMaxU = false;

    size_t heatSourceNumber = 2;
    size_t measrPointNumber = 4;
    double _lambda0 = +0.001;
    double _lambda1 = -0.01;
    double _factor = 1.0;

    DoubleMatrix alpha1;
    DoubleMatrix alpha2;
    DoubleMatrix alpha3;
    DoubleMatrix bettaX1;
    DoubleMatrix bettaX2;
    DoubleMatrix bettaX3;
    DoubleMatrix bettaY1;
    DoubleMatrix bettaY2;
    DoubleMatrix bettaY3;
    DoubleMatrix nomnU1;
    DoubleMatrix nomnU2;
    SpacePoint *measurePoints;

    ProblemParams *sourceParams = nullptr;

    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }

    DoubleMatrix V, U;

    HeatEquationIBVP forward;
    HeatEquationFBVP backward;

    void vectorToParameter(const DoubleVector &x);
    void parameterToVector(DoubleVector &x);

    bool drawImage = false;
    ODESolverMethod method = ODESolverMethod::EULER;
    GradientMethod *gradMethod;

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;

private:
    double _initialTemperature = NAN;
    double _environmentTemperature = NAN;

    double _initialTemperatureList[3] = { 0.5, 1.5, 2.0 };
    double _environmentTemperatureList[3] = { 0.5, 1.0, 1.5 };
    double _ratio1 = 1.0/1.0;
    double _ratio2 = 1.0/1.0;
    size_t _size1 = 1;
    size_t _size2 = 1;
};

}

#endif // PROBLEM3P_SOLVER1_H
