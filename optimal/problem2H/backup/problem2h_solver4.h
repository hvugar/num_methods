#ifndef PROBLEM2H_SOLVER4_H
#define PROBLEM2H_SOLVER4_H

#include "problem2h_common.h"
#include "problem2h_ibvp.h"

class PROBLEM2HSHARED_EXPORT Problem2HDirichlet4 : public RnFunction, public IGradient, public InitialBoundaryValueProblemPDE, public IProjection, public IPrinter, public IVectorNormalizer
{
public:
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HDirichlet4 &prob);
    static void checkGradient2(const Problem2HDirichlet4 &prob);

    Problem2HDirichlet4();
    virtual ~Problem2HDirichlet4();

    virtual auto fx(const DoubleVector &x) const -> double;
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;

    virtual auto norm(const DoubleVector &v) const -> double;
    virtual auto normalize(DoubleVector &v) const -> void;

    /** Integral part of functional */
    double integral(const std::vector<DoubleMatrix> &u) const;
    double integralU(const DoubleMatrix &u) const;
    double norm(const EquationParameterH &eprm, const OptimizeParameterH &oprm, const OptimizeParameterH &r_prm) const;

    /** Penalty part of functional */
    double penalty(const spif_vectorH &info, const OptimizeParameterH &o_prm) const;
    double gpi(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double g0i(unsigned int i, unsigned int layer, const spif_vectorH &u_info, const OptimizeParameterH &o_prm) const;
    double sign(double x) const;

public:
    void PrmToVector(const OptimizeParameterH &prm, DoubleVector &x) const;
    void VectorToPrm(const DoubleVector &x, OptimizeParameterH &prm) const;

    virtual double boundary(const SpaceNodePDE &, const TimeNodePDE &) const { return NAN; }

    virtual void print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const;

    virtual void project(DoubleVector &x, unsigned int index);
    virtual void project(DoubleVector &x) const;

    virtual void projectControlPoints(DoubleVector &x, unsigned int index) const;
    virtual void projectMeasurePoints(DoubleVector &x, unsigned int index) const;

    // ibvp -------------------------------------
    void initDeltaGrids(std::vector<DeltaGrid2D> &pulseDeltaGridList, std::vector<DeltaGrid2D> &msrntDeltaGridList, std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                        const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter) const;
    void releaseDeltaGrids(std::vector<DeltaGrid2D> &pulseDeltaGridList, std::vector<DeltaGrid2D> &msrntDeltaGridList, std::vector<DeltaGrid2D> &f_add2Info) const;
    void f_initialLayers(DoubleMatrix &u00, DoubleMatrix &u10, spif_vectorH &u_info, bool use,
                         unsigned int N, double hx, unsigned int M, double hy,
                         double ht, double aa__hxhx, double aa__hyhy, double lambda,
                         const std::vector<DeltaGrid2D> &pulseDeltaGridList,
                         const std::vector<DeltaGrid2D> &msrntDeltaGridList,
                         const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const;
    void b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vectorH &p_info, bool use,
                         unsigned int N, double hx, unsigned int M, double hy,
                         double ht, double aa__hxhx, double aa__hyhy, double lambda,
                         const std::vector<DeltaGrid2D> &pulseDeltaGridList,
                         const std::vector<DeltaGrid2D> &msrntDeltaGridList,
                         const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const;
    void f_borderCalculate(DoubleMatrix &u, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;
    void b_borderCalculate(DoubleMatrix &p, unsigned int N, double hx, unsigned int M, double hy, const TimeNodePDE &tn) const;
    void solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use) const;

    // ibvp -------------------------------------
    void solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info) const;
    void f_currentLayer(const DoubleMatrix &u10, const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter,
                        const std::vector<DeltaGrid2D> &msrntDeltaGridList, const std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                        unsigned int N, double hx, unsigned int M, double hy, DoubleMatrix &fxv, unsigned int ln) const;
    void b_currentLayer(const DoubleMatrix &p10, const EquationParameterH &equationParameter, const OptimizeParameterH &optimizeParameter,
                        const std::vector<DeltaGrid2D> &msrntDeltaGridList, const std::vector<DeltaGrid2D> &cntrlDeltaGridList,
                        unsigned int N, double hx, unsigned int M, double hy, DoubleMatrix &fxv, unsigned int ln, const spif_vectorH &u_info) const;
    //void b_initialLayers(DoubleMatrix &p00, DoubleMatrix &p10, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, unsigned int N, unsigned int M,
    //                     double hx, double hy, double ht, double aa__hxhx, double aa__hyhy, double lambda,
    //                     const std::vector<ExtendedSpacePointH> &cntExtSpacePoints,
    //                     const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const;
    double f_initial1(const SpaceNodePDE &sn) const;
    double b_initial1(const SpaceNodePDE &sn) const;
    double f_initial2(const SpaceNodePDE &sn) const;
    double b_initial2(const SpaceNodePDE &sn) const;
    double f_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;
    double b_boundary(const SpaceNodePDE &sn, const TimeNodePDE &tn) const;

    //void f_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
    //                    const std::vector<ExtendedSpacePointH> &cntExtSpacePoints, const std::vector<ExtendedSpacePointH> &msnExtSpacePoints) const;
    //void b_findRowsCols(GridH &grid, uint_vectorH &rows0, uint_vectorH &rows1, uint_vectorH &rows2, uint_vectorH &cols0, uint_vectorH &cols1, uint_vectorH &cols2, unsigned int N, unsigned int M,
    //                    const std::vector<ExtendedSpacePointH> &msnExtSpacePoints, const std::vector<ExtendedSpacePointH> &cntExtSpacePoints) const;
    auto prepareInfo(unsigned int N, const std::vector<SpacePoint> &points, spif_vectorH &info, unsigned int size) const -> void;
    void f_borderLayer(DoubleMatrix &u, DoubleMatrix &uh, unsigned int ln) const;
    auto f_add2Info(const DoubleMatrix &u, spif_vectorH &u_info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &msrntDeltaGridList) const -> void;
    auto b_add2Info(const DoubleMatrix &p, spif_vectorH &p_info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &cntrlDeltaGridList) const -> void;
    auto add2Info(const DoubleMatrix &u, spif_vectorH &info, unsigned int ln, double hx, double hy, const std::vector<DeltaGrid2D> &deltaList) const -> void;
    void f_layerInfo(const DoubleMatrix &u, unsigned int ln) const;
    void b_layerInfo(const DoubleMatrix &p, unsigned int ln) const;
    auto b_characteristic(const DoubleMatrix &u, unsigned int n, unsigned int m) const -> double;

    // common -----------------------------------
    //auto distributeTimeDelta(double t, double ht, unsigned int ln, const SpaceNodePDE &sn, const std::vector<ExtendedSpacePointH> &qPointNodes) const -> double;
    //auto newDistributeDeltaGaussPulse(const std::vector<SpacePoint> &thetas, std::vector<ExtendedSpacePointH> &extThetas, const Dimension &dimX, const Dimension &dimY) const -> void;
    //auto newDistributeDeltaGaussCntrl(const std::vector<SpacePoint> &cntrls, std::vector<ExtendedSpacePointH> &extCntrls, const Dimension &dimX, const Dimension &dimY) const -> void;
    //auto newDistributeDeltaGaussMsmnt(const std::vector<SpacePoint> &msmnts, std::vector<ExtendedSpacePointH> &extMsmnts, const Dimension &dimX, const Dimension &dimY) const -> void;

    auto v(unsigned int i, OptimizeParameterH o_prm, EquationParameterH e_prm, const spif_vectorH &u_info, unsigned int ln) const -> double;
    auto mu(unsigned int, unsigned int) const -> double { return 1.0; }
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

    GradientMethod *gm;
    std::vector<DoubleMatrix> vu;
    //DoubleMatrix f_initLayer;
    //DoubleMatrix b_initLayer;

    DoubleMatrix mPulseWeightMatrix;

    auto initiatePulseGrid() const -> void;
};

#endif // PROBLEM2H_SOLVER4_H
