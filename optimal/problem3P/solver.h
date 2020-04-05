#ifndef PROBLEM3P_SOLVER_H
#define PROBLEM3P_SOLVER_H

#include "global.h"

using namespace std;
using namespace std::placeholders;

namespace p3p
{

class Common;
class WaveEquationIBVP;
class WaveEquationFBVP;
class Solver;

struct ExternalSource
{
    SpacePoint point;
    double power;
    double sigmaX;
    double sigmaY;
    double sigmaT;
};

//struct FunctionalParemeter
//{
//    double epsilon1;
//    double epsilon2;
//};

struct Problem0HParameter
{
    SpacePoint p;
    std::vector<double> pwr_vl;
    std::vector<double> psi_vl;
    std::vector<double> psi_dx;
    std::vector<double> psi_dy;
    std::vector<double> psi_x;
    std::vector<double> psi_y;
    DeltaGrid2D deltaGrid;

    Problem0HParameter& initialize(const Dimension &time, const Dimension &dimX, const Dimension &dimY);
    void destroy();
    void distribute(const SpacePoint &p);
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP
{
public:
    HeatEquationIBVP(Solver *solver);
    virtual ~HeatEquationIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const override;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const override;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const override;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const override;
    virtual double weight() const override { return 0.5; }

public:
    virtual Dimension timeDimension() const override;
    virtual Dimension spaceDimensionX() const override;
    virtual Dimension spaceDimensionY() const override;
    virtual Dimension spaceDimensionZ() const override;

private:
    Solver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT FirstOrderLinearODE : virtual public IFirstOrderLinearODE
{
protected:
    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double;
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double;
    virtual auto B(const PointNodeODE &node, unsigned int row = 1) const -> double;
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const -> double;
    virtual auto boundary(const PointNodeODE &node, BoundaryConditionODE &condition, unsigned int row = 1) const -> double;
    virtual auto count() const -> unsigned int;

    virtual auto v(const PointNodeODE &node) const -> double;

    virtual auto z(const PointNodeODE &node, unsigned int row) const -> double;
    virtual auto zt(const PointNodeODE &node, unsigned int row) const -> double;
};

auto FirstOrderLinearODE::A(const PointNodeODE &, unsigned int row, unsigned int col) const -> double
{
    const double mx[2][2] = { {-3.0, -2.0}, { -4.0, -5.0 } };
    return mx[row-1][col-1];
}

auto FirstOrderLinearODE::B(const PointNodeODE &node, unsigned int row) const -> double
{
    return zt(node, row) - (A(node, row, 1)*z(node, 1)+A(node, row, 2)*z(node, 2)) - C(node, row)*v(node);
}

auto FirstOrderLinearODE::C(const PointNodeODE &, unsigned int row) const -> double
{
    const double c[2] = { +4.0, +1.0 };
    return c[row-1];
}

auto FirstOrderLinearODE::v(const PointNodeODE &node) const -> double
{
    return sin(2.0*M_PI*node.x*node.x);
}

auto FirstOrderLinearODE::initial(InitialCondition, unsigned int row) const -> double
{
    const double a[2] = { 0.4, 0.2 };
    return a[row-1];
}

auto FirstOrderLinearODE::boundary(const PointNodeODE &, BoundaryConditionODE &, unsigned int) const -> double
{
    throw new std::runtime_error("no boundary problem");
}

auto FirstOrderLinearODE::count() const -> unsigned int { return 2; }

auto FirstOrderLinearODE::z(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x;
    const double t2 = node.x*node.x;
    const double a[2] = { 0.5*sin(2.0*M_PI*t1)*sin(2.0*M_PI*t1) + 0.4*cos(2.0*M_PI*t2)*cos(2.0*M_PI*t2),
                          0.7*sin(4.0*M_PI*t1)*sin(4.0*M_PI*t1) + 0.2*cos(3.0*M_PI*t1)*cos(3.0*M_PI*t1)};
    return a[row-1];
}

auto FirstOrderLinearODE::zt(const PointNodeODE &node, unsigned int row) const -> double
{
    const double t1 = node.x;
    const double t2 = node.x*node.x;
    const double a[2] = { 0.5*M_PI*sin(2.0*M_PI*t1) - 1.6*M_PI*sin(4.0*M_PI*t2)*t1,
                          2.8*M_PI*sin(8.0*M_PI*t1) - 0.6*M_PI*sin(6.0*M_PI*t1)};
    return a[row-1];
}

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT HeatEquationFBVP : virtual public IHeatEquationFBVP
{
public:
    HeatEquationFBVP(Solver *solver);
    virtual ~HeatEquationFBVP();

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    virtual Dimension timeDimension() const;
    virtual Dimension spaceDimensionX() const;
    virtual Dimension spaceDimensionY() const;
    virtual Dimension spaceDimensionZ() const;

private:
    Solver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT FirstOrderLinearODE2 : virtual public IFirstOrderLinearODE
{
};

//--------------------------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT Solver : public RnFunction, public IGradient, public IProjection, public IPrinter, public R1Function
{
public:
    static void Main(int argc, char** argv);

    Solver();
    Solver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);
    virtual ~Solver();

private:
    Solver(const Solver &solver);

    static void optimization();

public:
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto fx(double t) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &) const  -> void;

    virtual auto print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f,
                       double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto integral1(const DoubleMatrix &u) const -> double;
    virtual auto integral2(const DoubleMatrix &u) const -> double;
    virtual auto norm() const -> double;
    virtual auto penalty() const -> double;

    void setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);

    auto vectorToParameter(const DoubleVector &x) const -> void;
    auto parameterToVector(DoubleVector &x) const -> void;

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    virtual void frw_calculate() const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    virtual void bcw_calculate() const;

    virtual void setTimeDimension(const Dimension &timeDimension) { this->_timeDimension = timeDimension; }
    virtual void setSpaceDimensionX(const Dimension &spaceDimensionX) { this->_spaceDimensionX = spaceDimensionX; }
    virtual void setSpaceDimensionY(const Dimension &spaceDimensionY) { this->_spaceDimensionY = spaceDimensionY; }

    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }
    virtual const Dimension& spaceDimensionZ() const { return _spaceDimensionZ; }

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;
    Dimension _spaceDimensionZ;

public:
    WaveEquationIBVP *forward = nullptr;
    WaveEquationFBVP *backward = nullptr;

    double p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    void frw_calculateU1U2(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void frw_saveToExcel(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void frw_saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void frw_saveToTextF(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    auto bcw_saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void;
    auto bcw_saveToImage(const DoubleMatrix &p, const TimeNodePDE &) const -> void;

    bool f_saveToFileTxt = false;
    bool f_saveToFilePng = false;

    ExternalSource external_source;
    double eps1 = 1.0;
    double eps2 = 1.0;

    DoubleMatrix u1;
    DoubleMatrix u2;

    unsigned int source_number = 0;

    std::vector<Problem0HParameter> optimalParameters;
};

};

#endif // PROBLEM0H_SOLVER_H
