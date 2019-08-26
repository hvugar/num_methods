#ifndef PROBLEM2H_SOLVER_H
#define PROBLEM2H_SOLVER_H

#include "problem2h_global.h"
#include <stdio.h>
#include <grid/hibvp.h>
#include <deltagrid.h>
#include <projection.h>
#include <gradient_sd.h>
#include <gradient_cjt.h>
#include <gradient_cs.h>

#define USE_NORM
#define USE_PENALTY

namespace Problem2H {

};

struct InitialPulse
{
    SpacePoint point;
    double blow;
};

struct PROBLEM2HSHARED_EXPORT SpacePointInfo
{
    SpacePointInfo(unsigned int length = 0);
    SpacePointInfo(const SpacePoint& point, unsigned int length);
    SpacePointInfo(const SpacePointInfo& spi);
    SpacePointInfo& operator= (const SpacePointInfo& other);
    virtual ~SpacePointInfo();

    void clear();

    SpacePoint point;
    unsigned int length;
    double *vl = nullptr;
    double *dx = nullptr;
    double *dy = nullptr;
    DeltaGrid2D deltaGrid;
};

//struct OptimizedParameter
//{
//    OptimizedParameter(unsigned int Nc, unsigned int No);
//    OptimizedParameter(const DoubleVector &x, unsigned int Nc, unsigned int No);

//    void frVector(const DoubleVector &x, unsigned int Nc, unsigned int No);
//    void toVector(DoubleVector &x, unsigned int Nc, unsigned int No);

//    const DoubleMatrix& k() const { return _k; }
//    DoubleMatrix& k() { return _k; }
//    const DoubleMatrix& z() const { return _z; }
//    DoubleMatrix& z() { return _z; }

//protected:
//    unsigned int Nc = 0;
//    unsigned int No = 0;
//    SpacePoint *eta = nullptr;
//    SpacePoint *ksi = nullptr;
//    DoubleMatrix _k;
//    DoubleMatrix _z;
//};

//struct RegulirizationParameter : public OptimizedParameter
//{
//    RegulirizationParameter(unsigned int Nc, unsigned int No);
//    RegulirizationParameter(const DoubleVector &x, unsigned int Nc, unsigned int No);

//    double regEpsilon1 = 0.0;
//    double regEpsilon2 = 0.0;
//    double regEpsilon3 = 0.0;
//    double regEpsilon4 = 0.0;
//};

struct Problem2HSharedData
{
protected:
    unsigned int Nc = 0;
    unsigned int No = 0;

    DoubleMatrix k;
    DoubleMatrix z;
    SpacePoint *eta = nullptr;
    SpacePoint *ksi = nullptr;

    DoubleMatrix r_k;
    DoubleMatrix r_z;
    SpacePoint *r_eta = nullptr;
    SpacePoint *r_ksi = nullptr;

public:
    double regEpsilon1 = 0.0;
    double regEpsilon2 = 0.0;
    double regEpsilon3 = 0.0;
    double regEpsilon4 = 0.0;

    bool optimizeK = true;
    bool optimizeZ = true;
    bool optimizeO = true;
    bool optimizeC = true;

    double *vmin = nullptr, *vmax = nullptr;
    double r = 1.0;

protected:
    SpacePointInfo *eta_info = nullptr;
    SpacePointInfo *ksi_info = nullptr;
    std::vector<DoubleMatrix> u_list;
};

class PROBLEM2HSHARED_EXPORT Problem2HCommon : virtual public Problem2HSharedData
{
public:
    virtual ~Problem2HCommon();

    virtual void setParameters(unsigned int Nc, unsigned int No, const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension);
    virtual void setOptimizationVector(const double *x);
    virtual void setRegularizationVector(const double *x, double espsilon1, double epsilon2, double epsilon3, double epsilon4);

    virtual auto penalty() const -> double;
    virtual auto gpi(unsigned int i, unsigned int ln) const -> double;
    virtual auto g0i(unsigned int i, unsigned int ln) const -> double;


    unsigned int Nt;
    TimeNodePDE *times;

    unsigned int L;
    unsigned int D;
};

class PROBLEM2HSHARED_EXPORT Problem2HWaveEquationIBVP :
        virtual public IWaveEquationIBVP,
        virtual public Problem2HCommon
{
protected:
    virtual double initial(const SpaceNodePDE &sn, InitialCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    virtual void setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY);
    void setInitialConditionMatrix(InitialPulse *initialPulses, unsigned int initialPulsesCount);
    void clrInitialConditionMatrix();

    unsigned int initialPulsesCount = 0;
    InitialPulse *initialPulses = nullptr;
    double noise = 0.0;
    bool save = false;
    std::vector<DoubleMatrix> save_u;

private:
    void layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);
    void layerInfoSave2TextFile(const DoubleMatrix &, const TimeNodePDE &) const;
    DoubleMatrix f_initialMatrix;
    DoubleMatrix f_crLayerMatrix;
    bool f_return_zero;
};

class PROBLEM2HSHARED_EXPORT Problem2HConjugateWaveEquationIBVP :
        virtual public IFinalWaveEquationIBVP,
        virtual public Problem2HCommon
{
protected:
    virtual double final(const SpaceNodePDE &sn, FinalCondition condition) const;
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn, BoundaryConditionPDE &condition) const;
    virtual double f(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    virtual void layerInfo(const DoubleMatrix &, const TimeNodePDE &) const;

public:
    virtual void setSpaceDimensions(const Dimension &dimensionX, const Dimension &dimensionY);

private:
    void layerInfoPrepareLayerMatrix(const DoubleMatrix &u, const TimeNodePDE& tn);
    void layerInfoSave2TextFile(const DoubleMatrix &u, const TimeNodePDE & tn) const;
    DoubleMatrix b_crLayerMatrix;
    bool f_return_zero;
};

class PROBLEM2HSHARED_EXPORT Problem2HSolver : virtual protected Problem2HWaveEquationIBVP,
        virtual protected Problem2HConjugateWaveEquationIBVP,
        virtual protected RnFunction,
        virtual protected IGradient,
        virtual protected IPrinter,
        virtual protected IProjection
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient3(const Problem2HSolver &prob, const double *data, unsigned int length);
    static void example1();
    static void example1_1(Problem2HSolver& ps, DoubleVector &x);
    static void example1_2(Problem2HSolver& ps, DoubleVector &x);

    static void example2();

    void setDimensions(const Dimension &dimensionX, const Dimension &dimensionY, const Dimension &timeDimension);
    void setEquationParameters(double waveSpeed, double waveDissipation);

    Problem2HWaveEquationIBVP& fw()   { return *(dynamic_cast<Problem2HWaveEquationIBVP*>(this)); }
    Problem2HConjugateWaveEquationIBVP& bw()  { return *(dynamic_cast<Problem2HConjugateWaveEquationIBVP*>(this)); }

protected:
    virtual double fx(const DoubleVector &x) const;
    virtual double fx_one(const DoubleVector &x, Problem2HSolver *solver) const;
    virtual double fx_norm() const;
    virtual void gradient(const DoubleVector &, DoubleVector &) const;
    virtual void gradient_one(const DoubleVector &, DoubleVector &, Problem2HSolver* solver) const;
    virtual void project(DoubleVector &x) const;
    virtual void project(DoubleVector &x, unsigned int);
    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;

    double integral(const std::vector<DoubleMatrix> &vu) const;
    double integralU(const DoubleMatrix &u) const;
};

#endif // PROBLEM2H_SOLVER_H
