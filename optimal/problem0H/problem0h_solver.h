#ifndef PROBLEM0H_SOLVER_H
#define PROBLEM0H_SOLVER_H

#include "problem0h_global.h"

using namespace std;
using namespace std::placeholders;

namespace h0p
{

class Problem0HCommon;
class WaveEquationIBVP;
class WaveEquationFBVP;
class ProblemSolver;

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

class PROBLEM0HSHARED_EXPORT WaveEquationIBVP : virtual public IWaveEquationIBVP
{
public:
    WaveEquationIBVP(ProblemSolver *solver);
    //WaveEquationIBVP(const WaveEquationIBVP &wibvp);
    virtual ~WaveEquationIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    virtual Dimension timeDimension() const;
    virtual Dimension spaceDimensionX() const;
    virtual Dimension spaceDimensionY() const;
    virtual Dimension spaceDimensionZ() const;

private:
    ProblemSolver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT WaveEquationFBVP : virtual public IWaveEquationFBVP
{
public:
    WaveEquationFBVP(ProblemSolver *solver);
    //WaveEquationFBVP(const WaveEquationFBVP &wibvp);
    virtual ~WaveEquationFBVP();

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
    ProblemSolver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT ProblemSolver : public RnFunction, public IGradient, public IProjection, public IPrinter, public R1Function
{
public:
    static void Main(int argc, char** argv);

    ProblemSolver();
    ProblemSolver(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);
    virtual ~ProblemSolver();

private:
    ProblemSolver(const ProblemSolver &solver);

    static void optimization();

public:
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto fx(double t) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    virtual auto project(DoubleVector &x, size_t index) -> void;
    virtual auto project(DoubleVector &) const  -> void;

    virtual auto print(unsigned int i, const DoubleVector &x, const DoubleVector &g, double f,
                       double alpha, GradientBasedMethod::MethodResult result) const -> void;

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
