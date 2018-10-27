#ifndef PROBLEM2P_SOLVER_H
#define PROBLEM2P_SOLVER_H

#include "problem2p_common.h"

class PROBLEM2PSHARED_EXPORT Problem2PNeumann : public RnFunction, public IGradient,
                                                public InitialBoundaryValueProblemPDE,
                                                public IProjection, public IPrinter,
                                                public IVectorNormalizer
{
public:
    static void Main(int argc, char** argv);
    static void checkGradient1(const Problem2PNeumann &prob, const OptimizeParameterP &o_prm);
    static void checkGradient2(const Problem2PNeumann &prob, const OptimizeParameterP &o_prm);

public:
    Problem2PNeumann();
    virtual ~Problem2PNeumann();

    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;
    virtual auto boundary(const SpaceNodePDE &, const TimeNodePDE &, BoundaryType) const -> double { return 0.0; }
    virtual auto project(DoubleVector &x, unsigned int index) -> void;
    virtual auto project(DoubleVector &x) const -> void;
    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g,
                       double f, double alpha, GradientMethod::MethodResult result) const -> void;
    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    /** Integral part of functional */
    auto mu(double x, double y) const -> double;
    auto fx_integral(const DoubleMatrix &u) const -> double;

    /** Norm part of functional */
    auto fx_norm(const EquationParameterP &eprm, const OptimizeParameterP &oprm, const OptimizeParameterP &rprm) const -> double;

    /** Penalty part of functional */
    auto fx_penalty(const spif_vector &info, const OptimizeParameterP &o_prm) const ->double;
    auto gpi(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameterP &o_prm) const -> double;
    auto g0i(unsigned int i, unsigned int layer, const spif_vector &u_info, const OptimizeParameterP &o_prm) const -> double;
    auto sign(double x) const -> double;

    auto projectControlPoints(DoubleVector &x, unsigned int index) const -> void;
    auto projectMeasurePoints(DoubleVector &x, unsigned int index) const -> void;

    /* Initial boundary value problems */
    auto solveForwardIBVP(DoubleMatrix &u, spif_vector &u_info, bool use, const OptimizeParameterP &mOptParameter) const -> void;
    auto solveBackwardIBVP(const DoubleMatrix &u, spif_vector &p_info, bool use, const spif_vector &u_info, const OptimizeParameterP &mOptParameter) const -> void;

    auto f_initialLayers(DoubleMatrix &u00, spif_vector &u_info, bool use, unsigned int N, unsigned int M,
                         double hx, double hy, const std::vector<ExtendedSpacePointP> &msnExtSpacePoints) const -> void;
    auto b_initialLayers(DoubleMatrix &p00, spif_vector &p_info, bool use, unsigned int N, unsigned int M,
                         double hx, double hy, const std::vector<ExtendedSpacePointP> &cntExtSpacePoints,
                         const DoubleMatrix &u) const -> void;
    auto f_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2,
                        uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                        unsigned int N, unsigned int M,
                        const std::vector<ExtendedSpacePointP> &cntExtSpacePoints,
                        const std::vector<ExtendedSpacePointP> &msnExtSpacePoints) const -> void;
    auto b_findRowsCols(uint_vector &rows0, uint_vector &rows1, uint_vector &rows2,
                        uint_vector &cols0, uint_vector &cols1, uint_vector &cols2,
                        unsigned int N, unsigned int M,
                        const std::vector<ExtendedSpacePointP> &msnExtSpacePoints,
                        const std::vector<ExtendedSpacePointP> &cntExtSpacePoints) const -> void;
    auto f_prepareInfo(unsigned int No, const std::vector<SpacePoint> &points, spif_vector &u_info, unsigned int L) const -> void;
    auto b_prepareInfo(unsigned int Nc, const std::vector<SpacePoint> &points, spif_vector &p_info, unsigned int L) const -> void;

    auto f_initial(const SpaceNodePDE &sn) const -> double;
    auto b_initial(const SpaceNodePDE &sn, const DoubleMatrix &u) const -> double;

    auto f_add2Info(const DoubleMatrix &u, spif_vector &u_info, unsigned int ln, double hx, double hy,
                    const std::vector<ExtendedSpacePointP> &extMsmnts, int method = 4) const -> void;
    auto b_add2Info(const DoubleMatrix &p, spif_vector &p_info, unsigned int ln, double hx, double hy,
                    const std::vector<ExtendedSpacePointP> &extCntrls, int method = 4) const -> void;

    auto f_layerInfo(const DoubleMatrix &u, unsigned int ln) const -> void;
    auto b_layerInfo(const DoubleMatrix &p, unsigned int ln) const -> void;

    auto newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePointP> &extCntrls, const Dimension &dimX, const Dimension &dimY) const -> void;
    auto newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePointP> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const -> void;

    auto VectorToPrm(const DoubleVector &pv, OptimizeParameterP &prm) const -> void;
    auto PrmToVector(const OptimizeParameterP &prm, DoubleVector &pv) const -> void;

    EquationParameterP mEquParameter;
    OptimizeParameterP mRegParameter;
    double r;

    DoubleVector vmin;
    DoubleVector vmax;

    double regEpsilon;
    DoubleMatrix U;

    bool optimizeK;
    bool optimizeZ;
    bool optimizeC;
    bool optimizeO;
};

#endif // PROBLEM2P_SOLVER_H
