#ifndef PROBLEM2H_SOLVER_BASE_H
#define PROBLEM2H_SOLVER_BASE_H

#include "problem2h_common.h"

class Problem2HDirichletBase : public InitialBoundaryValueProblemPDE, public RnFunction, public IGradient, public IProjection, public IPrinter, public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HDirichletBase &prob);
    static void checkGradient2(const Problem2HDirichletBase &prob);

    Problem2HDirichletBase();
    virtual ~Problem2HDirichletBase();

    /** Functional and gradient **/
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void = 0;
    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleMatrix> &u) const -> double = 0;
    virtual auto integralU(const DoubleMatrix &u) const -> double;
    /** Penalty part of functional */
    virtual auto penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const -> double = 0;
    virtual auto gpi(unsigned int i, unsigned int ln, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double;
    virtual auto g0i(unsigned int i, unsigned int ln, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const -> double;
    /** Norm part of functional */
    virtual auto norm(const EquationParameterH &eprm, const OptimizeParameterH &oprm, const OptimizeParameterH &r_prm) const -> double;
    /** ibv **/
    virtual auto solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &x, double lambda=0.25) const -> void = 0;
    virtual auto solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &x, double lambda=0.25) const -> void = 0;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto initPulseWeightMatrix(const std::vector<InitialPulse> &pulses) const -> void;

    virtual auto f_initial1(const SpaceNodePDE &sn) const -> double;
    virtual auto f_initial2(const SpaceNodePDE &sn) const -> double;
    virtual auto f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

    virtual auto b_initial1(const SpaceNodePDE &sn) const -> double;
    virtual auto b_initial2(const SpaceNodePDE &sn) const -> double;
    virtual auto b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

    virtual auto f_layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void;
    virtual auto b_layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &x) const -> void;

    virtual auto projectControlPoints(DoubleVector &x, unsigned int index) const -> void;
    virtual auto projectMeasurePoints(DoubleVector &x, unsigned int index) const -> void;

    auto prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void;
    auto add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void;

    auto b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double;

    auto PrmToVector(const OptimizeParameterH &prm, DoubleVector &x) const -> void;
    auto VectorToPrm(const DoubleVector &x, OptimizeParameterH &prm) const -> void;
    auto v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double;
    auto mu(unsigned int, unsigned int) const -> double { return 1.0; }
    auto sign(double x) const -> double;

    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    auto currentLayerFGrid(const DoubleMatrix &u, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids,
                           unsigned int ln) const -> void;
    auto currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids,
                           unsigned int ln, const spif_vectorH &u_info) const -> void;

    auto setGridDimensions(const Dimension &time, const Dimension &dimX, const Dimension &dimY) -> void;
public:
    EquationParameterH mEquParameter;
    OptimizeParameterH mOptParameter;
    OptimizeParameterH mRegParameter;

    DoubleVector vmin;
    DoubleVector vmax;

    unsigned int LD;

    bool optimizeK = true;
    bool optimizeZ = true;
    bool optimizeC = true;
    bool optimizeO = true;

    double r = 0.0;
    double regEpsilon = 0.0;
    double noise = 0.0;

    std::vector<DoubleMatrix> vu;

protected:
    virtual double boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

protected:
    DoubleMatrix mPulseWeightMatrix;
    DoubleMatrix mCrFfxWeightMatrix;

    DoubleMatrix mCrBfxWeightMatrix;
};

#endif // PROBLEM2H_SOLVER_BASE_H
