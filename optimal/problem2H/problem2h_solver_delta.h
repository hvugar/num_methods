#ifndef PROBLEM2H_SOLVER_DELTA_H
#define PROBLEM2H_SOLVER_DELTA_H

#include "problem2h_common.h"

class PROBLEM2HSHARED_EXPORT Problem2HDirichletDelta :
        public InitialBoundaryValueProblemPDE,
        public RnFunction, public IGradient, public IProjection,
        public IPrinter//, public IVectorNormalizer
{
public:
    /* static methods */
    static void Main(int argc, char* argv[]);
    static void checkGradient1(const Problem2HDirichletDelta &prob);
    static void checkGradient2(const Problem2HDirichletDelta &prob);
    static void example1();

    Problem2HDirichletDelta();
    virtual ~Problem2HDirichletDelta();

    /** Functional**/
    virtual auto fx(const DoubleVector &x) const -> double;
    /** Integral part of functional */
    virtual auto integral(const std::vector<DoubleMatrix> &u) const -> double;
    virtual auto integralU(const DoubleMatrix &u) const -> double;
    /** Penalty part of functional */
    virtual auto penalty(const spif_vectorH &info, const EquaParameter2H &o_prm) const -> double;
    virtual auto gpi(unsigned int i, unsigned int s, const spif_vectorH &u_info, const EquaParameter2H &equaPrm) const -> double;
    virtual auto g0i(unsigned int i, unsigned int s, const spif_vectorH &u_info, const EquaParameter2H &equaPrm) const -> double;
    /** Norm part of functional */
    virtual auto norm(const EquaParameter2H &equaPrm) const -> double;
    /** Functional gradient **/
    virtual auto gradient(const DoubleVector &, DoubleVector &) const -> void;
    /** IPrinter **/
    virtual auto print(unsigned int iteration, const DoubleVector &x, const DoubleVector &g, double f, double alpha, GradientMethod::MethodResult result) const -> void;

    auto solveForwardIBVP(std::vector<DoubleMatrix> &u, spif_vectorH &u_info, bool use, const DoubleVector &pv, double lambda=0.25) const -> void;
    auto solveBackwardIBVP(const std::vector<DoubleMatrix> &u, spif_vectorH &p_info, bool use, const spif_vectorH &u_info, const DoubleVector &pv, double lambda=0.25) const -> void;

    virtual auto initPulseWeightMatrix(const std::vector<InitialPulse2D> &pulses) const -> void;

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

    auto currentLayerFGrid(const DoubleMatrix &u, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln) const -> void;
    auto currentLayerBGrid(const DoubleMatrix &p, const std::vector<DeltaGrid2D> &controlDeltaGrids, const std::vector<DeltaGrid2D> &measurementDeltaGrids, unsigned int ln, const spif_vectorH &u_info) const -> void;

    virtual auto setGridDimensions(const Dimension &time, const Dimension &dimX, const Dimension &dimY) -> void;
    virtual auto v(unsigned int i, unsigned int s, const EquaParameter2H &equaPrm, const spif_vectorH &u_info) const -> double;
    virtual auto boundary(const SpaceNodePDE &, const TimeNodePDE &) const -> double;
    virtual auto mu(unsigned int, unsigned int) const -> double;
    virtual auto sign(double x) const -> double;

    EquaParameter2H equaPrm;
    FuncParameter2H funcPrm;

private:

    bool optimizeK = true;
    bool optimizeZ = true;
    bool optimizeC = true;
    bool optimizeO = true;

    unsigned int LD;
    std::vector<DoubleMatrix> vu;
    bool printLayers = false;
    double noise = 0.0;

protected:
    DoubleMatrix mPulseWeightMatrix;
    DoubleMatrix mCrFfxWeightMatrix;
    DoubleMatrix mCrBfxWeightMatrix;

};

#endif // PROBLEM2H_SOLVER1_H
