#ifndef PROBLEM2H_SOLVER1_H
#define PROBLEM2H_SOLVER1_H

#include "problem2h_common.h"

class PROBLEM2HSHARED_EXPORT Problem2HDirichlet1 : public RnFunction,
        public IGradient, public InitialBoundaryValueProblemPDE, public IProjection, public IPrinter,
        public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HDirichlet1 &prob);
    static void checkGradient2(const Problem2HDirichlet1 &prob);

    Problem2HDirichlet1();
    virtual ~Problem2HDirichlet1();

    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    /** Integral part of functional */
    double integral(const std::vector<DoubleMatrix> &u) const;
    double integralU(const DoubleMatrix &u) const;

    /** Norm part of functional */
    double norm(const EquationParameterH &eprm, const OptimizeParameterH &oprm, const OptimizeParameterH &r_prm) const;

    /** Penalty part of functional */
    double penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const;
    double gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double sign(double x) const;

public:
    virtual auto boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double { return NAN; }

    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void;

    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &x) const -> void;

    virtual auto projectControlPoints(DoubleVector &x, unsigned int index) const -> void;
    virtual auto projectMeasurePoints(DoubleVector &x, unsigned int index) const -> void;

    auto solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, double lambda=0.25) const -> void;
    auto solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const -> void;

    auto f_initial1(const SpaceNodePDE &sn) const -> double;
    auto b_initial1(const SpaceNodePDE &sn) const -> double;

    auto f_initial2(const SpaceNodePDE &sn) const -> double;
    auto b_initial2(const SpaceNodePDE &sn) const -> double;

    auto f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;
    auto b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const -> double;

    auto f_layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void;
    auto b_layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void;
    auto b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double;

    auto prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void;
    auto add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void;

    auto v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double;
    auto mu(unsigned int, unsigned int) const -> double { return 1.0; }
    auto PrmToVector(const OptimizeParameterH &prm, DoubleVector &x) const -> void;
    auto VectorToPrm(const DoubleVector &x, OptimizeParameterH &prm) const -> void;

public:
    EquationParameterH mEquParameter;
    OptimizeParameterH mOptParameter;
    OptimizeParameterH mRegParameter;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;
    double regEpsilon;
    double noise = 0.0;

    unsigned int LD;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;

    bool usePenalty;
    bool useNormal;

    std::vector<DoubleMatrix> vu;

    //GradientMethod *gm;

    DoubleMatrix mPulseWeightMatrix;
    DoubleMatrix mCrFfxWeightMatrix;
    DoubleMatrix mCrBfxWeightMatrix;

    auto initiatePulseGrid() const -> void;
    auto currentLayerFGrid(const DoubleMatrix &u, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids) const -> void;
    auto currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids,
                           double ln, const spif_vectorH &u_info) const -> void;
};

#endif // PROBLEM2H_SOLVER1_H
