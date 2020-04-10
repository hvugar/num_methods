#ifndef PROBLEM3P_SOLVER_H
#define PROBLEM3P_SOLVER_H

//ffmpeg -i %08d.png -framerate 160 -pattern_type sequence -vcodec mpeg4 -s 400x400 -start_number 0 -filter:v "setpts=0.5*PTS" test.avi

#include "global.h"

#define TIME_MIN 0
#define TIME_MAX 4000//1600//8000
#define DIMX_MAX 200//100//400
#define DIMY_MAX 200//100//400

#define TIME_STEP 0.0010//0.0025//0.0005
#define DIMX_STEP 0.005//0.01//0.0025
#define DIMY_STEP 0.005//0.01//0.0025

#define SIGMA_X 0.01
#define SIGMA_Y 0.01

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

class PROBLEM3P_SHARED_EXPORT HeatEquationIBVP : virtual public IHeatEquationIBVP, virtual public IFirstOrderLinearODE
{
public:
    HeatEquationIBVP(Solver *solver = nullptr);
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
    virtual auto initial(InitialCondition condition, unsigned int row = 1) const  -> double override;
    virtual auto count() const  -> unsigned int override;
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;
    virtual auto C(const PointNodeODE &node, unsigned int row = 1) const -> double;

public:
    Solver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT HeatEquationFBVP : virtual public IHeatEquationFBVP, virtual public IFirstOrderLinearODEFBVP
{
public:
    HeatEquationFBVP(Solver *solver);
    virtual ~HeatEquationFBVP();

    virtual auto final(const SpaceNodePDE &sn, FinalCondition condition) const -> double override;
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const -> double override;
    virtual auto f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double override;
    virtual auto layerInfo(const DoubleMatrix &, const TimeNodePDE &) const -> void override;
    virtual auto weight() const -> double override;
    virtual auto timeDimension() const -> Dimension override;
    virtual auto spaceDimensionX() const -> Dimension override;
    virtual auto spaceDimensionY() const -> Dimension override;
    virtual auto spaceDimensionZ() const -> Dimension override;

    virtual auto A(const PointNodeODE &node, unsigned int row = 1, unsigned int col = 1) const -> double override;
    virtual auto B(const PointNodeODE &node, unsigned int row = 1) const -> double override;
    virtual auto final(FinalCondition condition, unsigned int row = 1) const -> double override;
    virtual auto count() const  -> unsigned int override;
    virtual auto dimension() const -> Dimension override;
    virtual auto iterationInfo(const DoubleVector &v, const PointNodeODE &node) const -> void override;

private:
    Solver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM3P_SHARED_EXPORT Solver : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    static void Main(int argc, char** argv);
    static void example1();

    Solver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);

    Solver();
    virtual ~Solver();

private:
    Solver(const Solver &solver);

    static void optimization();

public:
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &) const  -> void;

    virtual auto print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f,
                       double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto integral(const DoubleMatrix &u) const -> double;
    virtual auto norm() const -> double;
    virtual auto penalty() const -> double;

    void setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);

    auto vectorToParameter(const DoubleVector &x) const -> void;
    auto parameterToVector(DoubleVector &x) const -> void;

    virtual double frw_initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double frw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double frw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void frw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void frw_calculate() const;

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    virtual void bcw_calculate() const;

    virtual void setTimeDimension(const Dimension &timeDimension);
    virtual void setSpaceDimensionX(const Dimension &spaceDimensionX);
    virtual void setSpaceDimensionY(const Dimension &spaceDimensionY);

    virtual const Dimension& timeDimension() const { return _timeDimension; }
    virtual const Dimension& spaceDimensionX() const { return _spaceDimensionX; }
    virtual const Dimension& spaceDimensionY() const { return _spaceDimensionY; }

    auto z(const PointNodeODE& node, unsigned int row) const -> double;
    auto zt(const PointNodeODE& node, unsigned int row) const -> double;
    auto v(const PointNodeODE &node) const -> double;

    void validate();

protected:
    Dimension _timeDimension;
    Dimension _spaceDimensionX;
    Dimension _spaceDimensionY;

public:
    HeatEquationIBVP *forward = nullptr;
    HeatEquationFBVP *backward = nullptr;

    double p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    void frw_calculateU1U2(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    void frw_saveToExcel(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    void frw_saveToImage(const DoubleMatrix &U, const TimeNodePDE &tn) const;
    void frw_saveToTextF(const DoubleMatrix &U, const TimeNodePDE &tn) const;

    auto bcw_saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void;
    auto bcw_saveToImage(const DoubleMatrix &p, const TimeNodePDE &) const -> void;

    bool f_saveToFileTxt = false;
    bool f_saveToFilePng = false;
    double epsilon = 1.0;

    std::vector<double> mq;
    DoubleVector mv;
    std::vector<SpacePoint> psi;
    std::vector<SpacePoint> mz;
    std::vector<SpacePoint> phi;

private:
    DoubleMatrix U;
    DoubleMatrix V;
    DeltaGrid2D deltaZ;
};

};

#endif // PROBLEM0H_SOLVER_H
