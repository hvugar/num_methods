#ifndef PROBLEM0H_SOLVER_H
#define PROBLEM0H_SOLVER_H

#include <float.h>
#include <limits.h>
#include <grid/hibvp.h>
#include <benchmark.h>
#include <function.h>
#include <deltagrid.h>
#include "problem0h_global.h"

#include <functional>

using namespace std;
using namespace std::placeholders;

namespace h0p
{

class Problem0HCommon;
class WaveEquationIBVP;
class WaveEquationFBVP;
class ProblemSolver;

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

    void create(const Dimension &time, const Dimension &dimX, const Dimension &dimY);
    void destroy();
    void distribute(const SpacePoint &p);
};

void Problem0HParameter::create(const Dimension &time, const Dimension &dimX, const Dimension &dimY)
{
    destroy();

    deltaGrid.initGrid(static_cast<unsigned int>(dimX.size()), dimX.step(), static_cast<unsigned int>(dimY.size()), dimY.step());

    const unsigned int L = 2*static_cast<unsigned int>(time.size());
    pwr_vl.resize(L+1, 0.0);
    psi_vl.resize(L+1, 0.0);
    psi_dx.resize(L+1, 0.0);
    psi_dy.resize(L+1, 0.0);
    psi_x.resize(L+1, 0.0);
    psi_y.resize(L+1, 0.0);

    //pwr_vl[0] = 0.0;
    //pwr_vl[1] = 0.0;
    //pwr_vl[L] = 0.0;
}

void Problem0HParameter::destroy()
{
    deltaGrid.cleanGrid();
    pwr_vl.clear();
    psi_vl.clear();
    psi_dx.clear();
    psi_dy.clear();
    psi_x.clear();
    psi_y.clear();
}

void Problem0HParameter::distribute(const SpacePoint &p)
{
    this->p = p;
    deltaGrid.distributeGauss(p, 1, 1);
}

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT WaveEquationIBVP : virtual public IWaveEquationIBVP
{
public:
    WaveEquationIBVP(ProblemSolver *solver);
    virtual ~WaveEquationIBVP();

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

    virtual const Dimension& timeDimension() const;
    virtual const Dimension& spaceDimensionX() const;
    virtual const Dimension& spaceDimensionY() const;
    virtual const Dimension& spaceDimensionZ() const;

private:
    ProblemSolver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT WaveEquationFBVP : virtual public IWaveEquationFBVP
{
public:
    WaveEquationFBVP(ProblemSolver *solver);
    virtual ~WaveEquationFBVP();

protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

    virtual const Dimension& timeDimension() const;
    virtual const Dimension& spaceDimensionX() const;
    virtual const Dimension& spaceDimensionY() const;
    virtual const Dimension& spaceDimensionZ() const;

private:
    ProblemSolver *solver;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT ProblemSolver : public RnFunction, public IGradient, public IProjection, public IPrinter
{
public:
    static void Main(int argc, char** argv);

    ProblemSolver();
    virtual ~ProblemSolver();

private:
    ProblemSolver(const ProblemSolver &solver);

    WaveEquationIBVP *forward;
    WaveEquationFBVP *backward;

    static void compareGradients();
    static void checkingForwardProblem();
    static void optimization();

public:
    virtual auto fx(const DoubleVector &x) const -> double;
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

    virtual double bcw_final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double bcw_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double bcw_f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void bcw_layerInfo(const DoubleMatrix &u, const TimeNodePDE &tn) const;

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

private:
    double p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    void calculateU1U2(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToExcel(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToTextF(const DoubleMatrix &u, const TimeNodePDE &tn) const;

    auto saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void;
    auto bcw_saveToImage(const DoubleMatrix &p, const TimeNodePDE &) const -> void;

    bool f_saveToFileTxt = false;
    bool f_saveToFilePng = false;

    double a = 1.0;
    double gamma = 0.0;

    double epsilon1 = 1.0;
    double epsilon2 = 0.0;

    //double alpha1 = 1.0;
    //double alpha2 = 1.0;
    //double alpha3 = 1.0;
    SpacePoint ksi;

    //DoubleMatrix U1;
    //DoubleMatrix U2;
    DoubleMatrix u1;
    DoubleMatrix u2;

    DoubleMatrix uL0;
    DoubleMatrix uL1;
    DoubleMatrix uL2;

    double p_sigmaX;
    double p_sigmaY;
    double p_sigmaT;

    unsigned int c_sigmaX;
    unsigned int c_sigmaY;

    unsigned int source_number;

    std::vector<Problem0HParameter> optimalParameters;
};

};

#endif // PROBLEM0H_SOLVER_H
