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

class Problem0HCommon;
class Problem0HForward;
class Problem0HBckward;
class Problem0HFunctional;

struct Problem0HParameter
{
    SpacePoint p;
    std::vector<double> pwr_vl;
    std::vector<double> psi_vl;
    std::vector<double> psi_dx;
    std::vector<double> psi_dy;
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
    pwr_vl.resize(L+1, 5.0);
    psi_vl.resize(L+1, 0.0);
    psi_dx.resize(L+1, 0.0);
    psi_dy.resize(L+1, 0.0);

    pwr_vl[0] = 0.0;
    pwr_vl[1] = 0.0;
    pwr_vl[L] = 0.0;
}

void Problem0HParameter::destroy()
{
    deltaGrid.cleanGrid();
    pwr_vl.clear();
    psi_vl.clear();
    psi_dx.clear();
    psi_dy.clear();
}

void Problem0HParameter::distribute(const SpacePoint &p)
{
    this->p = p;
    deltaGrid.distributeGauss(p);
}

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT Problem0HCommon
{
public:
    Problem0HCommon();
    virtual ~Problem0HCommon();

public:
    //inline auto virtual mu1(const SpaceNodePDE &) const -> double { return 1.0; }
    //inline auto virtual mu2(unsigned int, unsigned int) const -> double { return 1.0; }

    double a = 1.0;
    double gamma = 0.0;

    double epsilon1 = 1.0;
    double epsilon2 = 0.0;

    double alpha1 = 1.0;
    double alpha2 = 1.0;
    double alpha3 = 1.0;
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

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT Problem0HForward : public CdIHyperbolicIBVP, public virtual Problem0HCommon
{
public:
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

private:
    double p(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    void calculateU1U2(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToExcel(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToImage(const DoubleMatrix &u, const TimeNodePDE &tn) const;
    void saveToTextF(const DoubleMatrix &u, const TimeNodePDE &tn) const;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT Problem0HBckward : public ConjugateCdIHyperbolicIBVP, public virtual Problem0HCommon
{
public:
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;


private:
    auto saveBackwardInformarion(const DoubleMatrix &p, const TimeNodePDE &tn) const -> void;
    auto saveToImage(const DoubleMatrix &p, const TimeNodePDE &) const -> void;
};

//--------------------------------------------------------------------------------------------------------------//

class PROBLEM0HSHARED_EXPORT Problem0HFunctional : public RnFunction, public IGradient,
        protected virtual Problem0HForward, protected virtual Problem0HBckward
{
public:
    static void Main(int argc, char** argv);

    static void compareGradients();
    static void checkingForwardProblem();

public:
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &x, DoubleVector &g) const -> void;

    auto forward(DoubleMatrix &u, double a, double gamma) const -> void;
    auto backward(DoubleMatrix &u, double a, double gamma) const -> void;

    virtual auto integral1(const DoubleMatrix &u) const -> double;
    virtual auto integral2(const DoubleMatrix &u) const -> double;
    virtual auto norm() const -> double;
    virtual auto penalty() const -> double;

    void setDimension(const Dimension &timeDimension, const Dimension &spaceDimensionX, const Dimension &spaceDimensionY);

    auto vectorToParameter(const DoubleVector &x) const -> void;
    auto parameterToVector(DoubleVector &x) const -> void;
};


#endif // PROBLEM0H_SOLVER_H
