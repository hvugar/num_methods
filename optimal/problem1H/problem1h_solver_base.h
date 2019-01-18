#ifndef PROBLEM1H_SOLVER_BASE_H
#define PROBLEM1H_SOLVER_BASE_H

#include "problem1h_common.h"

class PROBLEM1HSHARED_EXPORT Problem1HDirichletBase : public InitialBoundaryValueProblemPDE, public RnFunction, public IGradient, public IProjection, public IPrinter, public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem1HDirichletBase &prob);
    static void checkGradient2(const Problem1HDirichletBase &prob);

    Problem1HDirichletBase();
    virtual ~Problem1HDirichletBase();

    /** Functional and gradient **/
    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void = 0;
    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleVector> &u) const -> double = 0;
    virtual auto integralU(const DoubleVector &u) const -> double;
    /** Penalty part of functional */
    virtual auto penalty(const spif_vector1H &info, const OptimizeParameter1H &o_prm) const -> double = 0;
    virtual auto gpi(unsigned int i, unsigned int ln, const spif_vector1H &u_info, const OptimizeParameter1H &o_prm) const -> double;
    virtual auto g0i(unsigned int i, unsigned int ln, const spif_vector1H &u_info, const OptimizeParameter1H &o_prm) const -> double;
    /** Norm part of functional */
    virtual auto norm(const EquationParameter1H &eprm, const OptimizeParameter1H &oprm, const OptimizeParameter1H &r_prm) const -> double;
    /** ibv **/
    virtual auto solveForwardIBVP(std::vector<DoubleVector> &u, spif_vector1H &u_info, bool use, double lambda=0.25) const -> void = 0;
    virtual auto solveBackwardIBVP(const std::vector<DoubleVector> &u, spif_vector1H &p_info, bool use, const spif_vector1H &u_info, double lambda=0.25) const -> void = 0;

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto initPulseWeightVector(const std::vector<PulseSpacePoint1H> &theta) const -> void;

    virtual auto f_initial1(const SpaceNodePDE &sn) const -> double;
    virtual auto f_initial2(const SpaceNodePDE &sn) const -> double;
    virtual auto f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

    virtual auto b_initial1(const SpaceNodePDE &sn) const -> double;
    virtual auto b_initial2(const SpaceNodePDE &sn) const -> double;
    virtual auto b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

    virtual auto f_layerInfo(const DoubleVector &u, unsigned int ln) const -> void;
    virtual auto b_layerInfo(const DoubleVector &p, unsigned int ln) const -> void;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &x) const -> void;

    virtual auto projectControlPoints(DoubleVector &x, unsigned int index) const -> void;
    virtual auto projectMeasurePoints(DoubleVector &x, unsigned int index) const -> void;

    auto prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vector1H &info, unsigned int size) const -> void;
    auto add2Info(const DoubleVector &u, spif_vector1H &info, unsigned int ln, double hx, const std::vector<DeltaGrid1D> &deltaList) const -> void;

    auto b_characteristic(const DoubleVector &u, unsigned int n) const -> double;

    auto PrmToVector(const OptimizeParameter1H &prm, DoubleVector &x) const -> void;
    auto VectorToPrm(const DoubleVector &x, OptimizeParameter1H &prm) const -> void;
    auto v(unsigned int i, OptimizeParameter1H o_prm, EquationParameter1H e_prm, const spif_vector1H &u_info, unsigned int ln) const -> double;
    auto mu(unsigned int) const -> double { return 1.0; }
    auto sign(double x) const -> double;

    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    auto currentLayerFGrid(const DoubleVector &u, const std::vector<DeltaGrid1D> &controlDeltaGrids, const std::vector<DeltaGrid1D> &measurementDeltaGrids) const -> void;
    auto currentLayerBGrid(const DoubleVector &p, const std::vector<DeltaGrid1D> &controlDeltaGrids, const std::vector<DeltaGrid1D> &measurementDeltaGrids,
                           unsigned int ln, const spif_vector1H &u_info) const -> void;

    auto setGridDimensions(const Dimension &time, const Dimension &dimX) -> void;
public:
    EquationParameter1H mEquParameter;
    OptimizeParameter1H mOptParameter;
    OptimizeParameter1H mRegParameter;

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

    std::vector<DoubleVector> vu;

    bool printLayers = false;

protected:
    virtual auto boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

protected:
    DoubleVector mPulseWeightVector;
    DoubleVector mCrFfxWeightMatrix;
    DoubleVector mCrBfxWeightMatrix;
};

#endif // PROBLEM1H_SOLVER_BASE_H
